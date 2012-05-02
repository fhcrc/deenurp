"""
Generate reference clusters for use in refpkg building
"""
import argparse
import csv
import logging
import shutil

from Bio import SeqIO
from romperroom import uclust

from .. import tax, wrap

def build_parser(p):
    p.add_argument('named_sequence_file', help="""Named sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""", type=argparse.FileType('r'))
    p.add_argument('taxtable', help="""Taxtable""", type=argparse.FileType('r'))
    p.add_argument('sequence_out', type=argparse.FileType('w'))
    p.add_argument('seqinfo_out', type=argparse.FileType('w'))
    p.add_argument('cluster_info_out', type=argparse.FileType('w'))
    p.add_argument('--cluster-rank', help="""Rank to cluster sequences
            [default: %(default)s]""", default='species')
    p.add_argument('-u', '--unnamed-sequences', help="""Path to unnamed sequence file""")
    p.add_argument('-r', '--redundant-cluster-id', default=0.995, type=float)
    p.add_argument('-i', '--cluster-id', default=0.985, type=float)

def cluster_identify_redundant(named_sequence_file, named_ids, to_cluster,
        threshold=0.995):
    with wrap.ntf(suffix='.uc', prefix='to_cluster') as tf:
        # Search with uclust
        uclust.search(named_sequence_file, to_cluster, tf.name,
                pct_id=threshold, search_pct_id=0.90, trunclabels=True)
        records = uclust.parse_uclust_out(tf)
        hits = (i.query_label for i in records if i.type == 'H')
        return frozenset(hits)

def cluster_info_writer(fp):
    w = csv.writer(fp, quoting=csv.QUOTE_NONNUMERIC, lineterminator='\n')
    w.writerow(('seqname', 'cluster'))
    def write_cluster_info(sequence_id, cluster_id):
        w.writerow((sequence_id, cluster_id))
    return write_cluster_info

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
        taxonomy = tax.TaxNode.from_taxtable(fp)
    with a.seqinfo_file as fp:
        logging.info('Loading seqinfo')
        r = csv.DictReader(fp)
        seqinfo = {i['seqname']: i for i in r}

    # Add sequences to taxonomy
    for i in seqinfo.values():
        taxonomy.get_node(i['tax_id']).sequence_ids.append(i['seqname'])

    # Write clustering information for sequences with cluster_rank-level
    # classifications
    done = set()
    with a.sequence_out, a.cluster_info_out as cinfo_out:
        w = cluster_info_writer(cinfo_out)
        for tax_id, sequences in taxonomic_clustered(taxonomy, a.cluster_rank):
            done |= sequences
            for s in sequences:
                w(s, tax_id)

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
        with wrap.ntf(prefix='to_cluster', suffix='.fasta') as tf, \
                wrap.ntf(prefix='unnamed_to_cluster', suffix='.fasta') as unnamed_fp:
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

            # Cluster remaining sequences into OTUs
            for i, cluster_seqs in enumerate(identify_otus_unnamed(unnamed_fp.name, a.cluster_id)):
                done |= cluster_seqs
                otu = 'otu_{0}'.format(i)
                for sequence in cluster_seqs:
                    w(sequence, otu)

    with a.seqinfo_out as fp:
        seqinfo_records = (seqinfo.get(i, {'seqname': i}) for i in done)
        w = csv.DictWriter(fp, seqinfo.values()[0].keys(),
                quoting=csv.QUOTE_NONNUMERIC, lineterminator='\n')
        w.writeheader()
        w.writerows(seqinfo_records)
