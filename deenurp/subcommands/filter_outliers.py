"""
Filter sequences at the species level, whose distance from the species medoid
is above some threshold.
"""

import argparse
import csv
import logging
import os.path
import subprocess

from Bio import SeqIO

from .. import wrap, tax

DEFAULT_RANK = 'species'
RSCRIPT_PATH = os.path.join(os.path.dirname(__file__),
        '..', 'data', 'find_outliers.R')

def build_parser(p):
    p.add_argument('sequence_file', help="""All sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""",
            type=argparse.FileType('r'))
    p.add_argument('taxonomy', help="""Taxtable""", type=argparse.FileType('r'))
    p.add_argument('output_fp', help="""Destination for sequences""",
            type=argparse.FileType('w'))
    p.add_argument('--filter-rank', default=DEFAULT_RANK)
    p.add_argument('--distance-cutoff', type=float, default=0.015,
            help="""Distance cutoff from cluster centroid [default:
            %(default)f]""")
    p.add_argument('--threads', type=int, default=12)

def load_seqinfo_into_tax(seqinfo_fp, taxonomy):
    r = csv.DictReader(seqinfo_fp)
    c = 0
    for i in r:
        taxonomy.get_node(i['tax_id']).sequence_ids.append(i['seqname'])
        c += 1
    return c

def sequences_above_rank(taxonomy, rank=DEFAULT_RANK):
    """
    Generate the sequence ids from taxonomy whose rank is above specified rank.
    """
    ranks = taxonomy.ranks
    r_index = ranks.index(rank)
    assert r_index >= 0
    def above_rank(node):
        n_index = ranks.index(node.rank)
        assert n_index >= 0
        return n_index < r_index

    for n in taxonomy:
        if above_rank(n):
            for sequence_id in n.sequence_ids:
                yield sequence_id

def run_r_find_outliers(sequence_file, cutoff):
    """
    Run the R script find_outliers.R
    """
    with wrap.ntf(prefix='prune-') as tf:
        cmd = [RSCRIPT_PATH, sequence_file, str(cutoff), tf.name]
        subprocess.check_call(cmd)
        return [i.strip() for i in tf]

def filter_sequences(sequence_file, cutoff, threads=12):
    with wrap.ntf(prefix='cmalign', suffix='.sto') as a_sto, \
         wrap.ntf(prefix='cmalign', suffix='.fasta') as a_fasta, \
         open(os.devnull) as devnull:
        # Align
        wrap.cmalign_files(sequence_file, a_sto.name,
                stdout=devnull, mpi_args=['-np', str(threads)])
        # APE requires FASTA
        SeqIO.convert(a_sto, 'stockholm', a_fasta, 'fasta')
        a_fasta.flush()
        return run_r_find_outliers(a_fasta.name, cutoff)

def action(a):
    # Load taxonomy
    with a.taxonomy as fp:
        taxonomy = tax.TaxNode.from_taxtable(fp)
        logging.info('Loaded taxonomy')

    # Load sequences into taxonomy
    with a.seqinfo_file as fp:
        count = load_seqinfo_into_tax(fp, taxonomy)
        logging.info('Added %d sequences', count)

    # Sequences which are classified above the desired rank should just be kept
    kept_ids = frozenset(sequences_above_rank(taxonomy, a.filter_rank))

    with a.output_fp as fp:
        logging.info('Keeping %d sequences classified above %s', len(kept_ids), a.filter_rank)
        wrap.esl_sfetch(a.sequence_file, kept_ids, fp)

        # For each filter-rank, filter
        nodes = [i for i in taxonomy if i.rank == a.filter_rank]
        keep_ids = set()
        for i, node in enumerate(nodes):
            seqs = frozenset(node.subtree_sequence_ids())
            if not seqs:
                logging.warn("No sequences for %s (%s)", node.tax_id, node.name)
                continue
            elif len(seqs) < 3:
                logging.warn("Only %d sequences for %s (%s). Keeping.", len(seqs),
                        node.tax_id, node.name)
                keep_ids |= frozenset(seqs)
                continue
            with wrap.ntf(prefix='to_filter', suffix='.fasta') as tf:
                logging.warn("%d sequences for %s (%s). [%d/%d]", len(seqs),
                        node.tax_id, node.name, i + 1, len(nodes))
                # Extract sequences
                wrap.esl_sfetch(a.sequence_file, seqs,
                        tf)
                tf.flush()
                prune = frozenset(filter_sequences(tf.name, a.distance_cutoff))
                logging.info("Pruning %d", len(prune))
                assert not prune - seqs
                keep_ids |= seqs - prune

                # Extract
                wrap.esl_sfetch(a.sequence_file, seqs - prune, fp)

        # Extract all hits
        #wrap.esl_sfetch(a.sequence_file, kept_ids, fp)
