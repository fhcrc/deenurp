"""
Generate reference clusters for use in refpkg building
"""
import argparse
import csv
import logging
import shutil

from Bio import SeqIO
from deenurp import uclust
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
    p.add_argument('-m', '--unnamed-sequence-meta', help="""Path to unnamed
            sequence file sequence info""", type=argparse.FileType('r'))
    p.add_argument('-r', '--redundant-cluster-id', default=0.97, type=float,
            help="""Percent ID at which to remove redundant sequences.
            [default: %(default).3f]""")
    p.add_argument('-i', '--cluster-id', default=0.985, type=float, help="""Cluster ID [default: %(default).3f]""")

def cluster_identify_redundant(named_sequence_file, named_ids, to_cluster,
        threshold=0.97):
    with util.ntf(suffix='.uc', prefix='to_cluster') as tf:
        # Search with uclust
        uclust.search(named_sequence_file, to_cluster, tf.name,
                pct_id=0.80,
                maxaccepts=5,
                maxrejects=100)

        # Uclust.search renames to tf, need a new handle.
        records = uclust.parse_uclust_out(tf)
        hits = (i.query_label for i in records
                if i.type == 'H' and i.pct_id >= threshold * 100.0)

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
    logging.info('Running UCLUST on unnamed sequences at %f',
                 cluster_similarity)
    with util.ntf(prefix='uclust') as tf:
        # Sort and cluster
        uclust.cluster(
            seq_file, tf.name, pct_id=cluster_similarity, quiet=True)
        clusters = uclust.sequences_by_cluster(uclust.parse_uclust_out(tf))
        for _, sequences in clusters:
            yield [i.query_label for i in sequences]

def action(a):
    # index fasta file
    fa_idx = wrap.read_seq_file(a.named_sequence_file)

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
    if a.unnamed_sequence_meta:
        with a.unnamed_sequence_meta as fp:
            r = csv.DictReader(fp)
            unnamed_seqinfo = {i['seqname']: i for i in r}
        assert not set(unnamed_seqinfo) & set(seqinfo)
        seqinfo.update(unnamed_seqinfo)

    # Write clustering information for sequences with cluster_rank-level
    # classifications
    done = set()
    cluster_ids = {}
    with a.sequence_out:
        for tax_id, sequences in taxonomic_clustered(taxonomy, a.cluster_rank):
            for sequence in sequences:
                cluster_ids[sequence] = tax_id
            done |= set(sequences)

        # Fetch sequences
        logging.info('Fetching %d %s-level sequences', len(done), a.cluster_rank)
        wrap.esl_sfetch(a.named_sequence_file, done, a.sequence_out, fa_idx)
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
            wrap.esl_sfetch(a.named_sequence_file, above_rank_seqs, tf, fa_idx)
            if a.unnamed_sequences:
                with open(a.unnamed_sequences) as fp:
                    shutil.copyfileobj(fp, tf)
            tf.close()

            # Remove redundant sequences: we don't need anything that's unnamed
            # & close to something named.
            redundant_ids = cluster_identify_redundant(a.sequence_out.name,
                    done, to_cluster=tf.name, threshold=a.redundant_cluster_id)
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
                done |= set(cluster_seqs)
                otu = 'otu_{0}'.format(i)
                for sequence in cluster_seqs:
                    cluster_ids[sequence] = otu

    with a.seqinfo_out as fp:
        def add_cluster(i):
            """Add a cluster identifier to sequence metadata"""
            i['cluster'] = cluster_ids[i['seqname']]
            return i
        seqinfo_records = (seqinfo.get(i, {'seqname': i}) for i in done)
        seqinfo_records = (add_cluster(i) for i in seqinfo_records)

        fields = list(list(seqinfo.values())[0].keys())
        fields.append('cluster')
        w = csv.DictWriter(fp, fields,
                quoting=csv.QUOTE_NONNUMERIC, lineterminator='\n')
        w.writeheader()
        w.writerows(seqinfo_records)
