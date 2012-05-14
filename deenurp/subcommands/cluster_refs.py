"""
Generate reference clusters for use in refpkg building
"""
import argparse
import csv
import logging
import shutil

from Bio import SeqIO
from romperroom import uclust
from taxtastic.taxtable import TaxNode

from .. import wrap, util

def build_parser(p):
    p.add_argument('named_sequence_file', help="""Named sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""", type=argparse.FileType('r'))
    p.add_argument('taxtable', help="""Taxtable""", type=argparse.FileType('r'))
    p.add_argument('sequence_out', type=argparse.FileType('w'))
    p.add_argument('seqinfo_out', type=argparse.FileType('w'))
    p.add_argument('--cluster-rank', help="""Rank to cluster sequences
            [default: %(default)s]""", default='species')
    p.add_argument('-u', '--unnamed-sequences', help="""Path to unnamed sequence file""")
    p.add_argument('-r', '--redundant-cluster-id', default=0.995, type=float)
    p.add_argument('-i', '--cluster-id', default=0.985, type=float)

def cluster_identify_redundant(named_sequence_file, named_ids, to_cluster,
        threshold=0.995):
    with util.ntf(suffix='.uc', prefix='to_cluster') as tf:
        # Search with uclust
        uclust.search(named_sequence_file, to_cluster, tf.name,
                pct_id=threshold, search_pct_id=0.90, trunclabels=True)
        records = uclust.parse_uclust_out(tf)
        hits = (i.query_label for i in records if i.type == 'H')
        return frozenset(hits)

def taxonomic_clustered(taxonomy, cluster_rank):
    """
    Generate tax_id, sequence_id_set tuples for each tax_id at cluster_rank
    """
    nodes = (node for node in taxonomy if node.rank == cluster_rank)
    return ((node.tax_id, frozenset(node.subtree_sequence_ids()))
            for node in nodes)

def identify_otus_unnamed(seq_file, cluster_similarity):
    """
    Generates sequence ids in a cluster

    Identify sequences in OTUs at the given cluster similarity;
    """
    for cluster in wrap.dnaclust(seq_file, similarity=cluster_similarity):
        yield cluster.sequences

def action(a):
    # Load taxtable
    with a.taxtable as fp:
        logging.info('Loading taxonomy')
        taxonomy = TaxNode.from_taxtable(fp)
    with a.seqinfo_file as fp:
        logging.info('Loading seqinfo')
        taxonomy.populate_from_seqinfo(fp)
        fp.seek(0)
        r = csv.DictReader(fp)
        seqinfo = {i['seqname']: i for i in r}

    # Write clustering information for sequences with cluster_rank-level
    # classifications
    done = set()
    with a.sequence_out:
        for tax_id, sequences in taxonomic_clustered(taxonomy, a.cluster_rank):
            done |= sequences

        # Fetch sequences
        logging.info('Fetching %d %s-level sequences', len(done), a.cluster_rank)
        wrap.esl_sfetch(a.named_sequence_file, done, a.sequence_out)
        a.sequence_out.flush()

        # Find sequences *above* cluster_rank
        above_rank_seqs = frozenset(i for i in taxonomy.subtree_sequence_ids()
                                    if i not in done)
        logging.info('%d sequences above rank %s', len(above_rank_seqs),
                a.cluster_rank)

        # Write sequences clustered above species level, unnamed sequences to
        # file
        with util.ntf(prefix='to_cluster', suffix='.fasta') as tf, \
                util.ntf(prefix='unnamed_to_cluster', suffix='.fasta') as unnamed_fp:
            wrap.esl_sfetch(a.named_sequence_file, above_rank_seqs, tf)
            if a.unnamed_sequences:
                with open(a.unnamed_sequences) as fp:
                    shutil.copyfileobj(fp, tf)
            tf.close()

            # Remove redundant sequences: we don't need anything that's unnamed
            # & close to something named.
            redundant_ids = cluster_identify_redundant(a.sequence_out.name,
                    done, tf.name, a.redundant_cluster_id)
            logging.info('%d redundant sequences', len(redundant_ids))

            # Extract desired sequences
            sequences = SeqIO.parse(tf.name, 'fasta')
            sequences = (i for i in sequences if i.id not in redundant_ids)

            # Write to file for clustering
            unnamed_count = SeqIO.write(sequences, unnamed_fp, 'fasta')
            logging.info('Kept %d non-redundant, unnamed sequences', unnamed_count)

            # Write to output sequence file
            unnamed_fp.seek(0)
            shutil.copyfileobj(unnamed_fp, a.sequence_out)
            unnamed_fp.close()

            cluster_ids = {}
            # Cluster remaining sequences into OTUs
            for i, cluster_seqs in enumerate(identify_otus_unnamed(unnamed_fp.name, a.cluster_id)):
                done |= cluster_seqs
                otu = 'otu_{0}'.format(i)
                for sequence in cluster_seqs:
                    cluster_ids[sequence] = otu

    with a.seqinfo_out as fp:
        def add_cluster(i):
            """Add a cluster identifier to sequence metadata"""
            if 'tax_id' in i:
                i['cluster'] = i['tax_id']
            else:
                i['cluster'] = cluster_ids[i['seqname']]
            return i
        seqinfo_records = (seqinfo.get(i, {'seqname': i}) for i in done)
        seqinfo_records = (add_cluster(i) for i in seqinfo_records)

        fields = list(seqinfo.values()[0].keys())
        fields.append('cluster')
        w = csv.DictWriter(fp, fields,
                quoting=csv.QUOTE_NONNUMERIC, lineterminator='\n')
        w.writeheader()
        w.writerows(seqinfo_records)
