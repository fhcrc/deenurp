"""
Filter sequences at the species level, whose distance from the species medoid
is above some threshold.
"""

import argparse
import csv
import logging
import os.path
import sys

from concurrent import futures

from Bio import SeqIO
from taxtastic.taxtable import TaxNode

from .. import wrap, util, outliers

DEFAULT_RANK = 'species'
RSCRIPT_PATH = os.path.join(os.path.dirname(__file__),
        '..', 'data', 'find_outliers.R')

DROP = 'drop'
KEEP = 'keep'

def build_parser(p):
    p.add_argument('sequence_file', help="""All sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""",
            type=argparse.FileType('r'))
    p.add_argument('taxonomy', help="""Taxtable""", type=argparse.FileType('r'))
    p.add_argument('output_fp', help="""Destination for sequences""",
            type=argparse.FileType('w'))
    p.add_argument('--filter-rank', default=DEFAULT_RANK)
    p.add_argument('--filtered-seqinfo',
            help="""Path to write filtered sequence info""",
            type=argparse.FileType('w'))
    p.add_argument('--log', help="""Log path""", type=argparse.FileType('w'))
    p.add_argument('--distance-cutoff', type=float, default=0.015,
            help="""Distance cutoff from cluster centroid [default:
            %(default)f]""")
    p.add_argument('--threads', type=int, default=12)

    rare_group = p.add_argument_group("Rare taxa")
    rare_group.add_argument('--min-seqs-for-filtering', type=int, default=3, help="""Minimum
            number of sequences perform distance-based medoid-filtering on [default:
            %(default)d]""")
    rare_group.add_argument('--rare-taxon-action', choices=(KEEP, DROP), default=KEEP,
            help="""Action to perform when a taxon has < '--min-seqs-to-filter'
            representatives. [default: %(default)s]""")

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

def filter_sequences(sequence_file, tax_id, cutoff):
    """
    Return a list of sequence names identifying outliers.
    """

    prefix = '{}_'.format(tax_id)

    with util.ntf(prefix=prefix, suffix='.sto') as a_sto, \
            util.ntf(prefix=prefix, suffix='.fasta') as a_fasta:
        # Align
        wrap.cmalign_files(sequence_file, a_sto.name, threads=1)
        # FastTree requires FASTA
        SeqIO.convert(a_sto, 'stockholm', a_fasta, 'fasta')
        a_fasta.flush()

        taxa, distmat = outliers.fasttree_dists(a_fasta.name)
        is_out = outliers.outliers(distmat, cutoff)

        assert len(is_out) == len(taxa)

        return [t for t, o in zip(taxa, is_out) if o]

def filter_worker(sequence_file, node, seqs, distance_cutoff, log_taxid=None):
    """
    Worker task for running filtering tasks.

    Arguments:
    :sequence_file: Complete sequence file
    :node: Taxonomic node being filtered
    :seqs: set containing the names of sequences to keep
    :distance_cutoff: Distance cutoff for medoid filtering
    :sfetch_lock: Lock to acquire to retrieving sequences from ``sequence_file``
    :log_taxid: Optional function to log tax_id activity.

    :returns: Set of sequences to *keep*
    """

    prefix = '{}_'.format(node.tax_id)

    with util.ntf(prefix=prefix, suffix='.fasta') as tf:
        # Extract sequences
        wrap.esl_sfetch(sequence_file, seqs, tf)
        tf.flush()
        prune = frozenset(filter_sequences(tf.name, node.tax_id, distance_cutoff))
        assert not prune - seqs
        if log_taxid:
            log_taxid(node.tax_id, node.name, len(seqs), len(seqs - prune),
                      len(prune))
        return seqs - prune

def action(a):
    # Load taxonomy
    with a.taxonomy as fp:
        taxonomy = TaxNode.from_taxtable(fp)
        logging.info('Loaded taxonomy')

    # Load sequences into taxonomy
    with a.seqinfo_file as fp:
        taxonomy.populate_from_seqinfo(fp)
        logging.info('Added %d sequences', sum(1 for i in taxonomy.subtree_sequence_ids()))

    # Sequences which are classified above the desired rank should just be kept
    kept_ids = frozenset(sequences_above_rank(taxonomy, a.filter_rank))

    log_taxid = None
    if a.log:
        writer = csv.writer(a.log, lineterminator='\n',
                quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(('tax_id', 'tax_name', 'n', 'kept', 'pruned'))

        def log_taxid(tax_id, tax_name, n, kept, pruned):
            writer.writerow((tax_id, tax_name, n, kept, pruned))

    with a.output_fp as fp, a.log or util.nothing():
        logging.info('Keeping %d sequences classified above %s', len(kept_ids), a.filter_rank)

        # For each filter-rank, filter
        nodes = [i for i in taxonomy if i.rank == a.filter_rank]

        # Filter each tax_id, running in ``--threads`` tasks in parallel
        with futures.ThreadPoolExecutor(a.threads) as executor:
            futs = {}
            for i, node in enumerate(nodes):
                seqs = frozenset(node.subtree_sequence_ids())
                if not seqs:
                    logging.debug("No sequences for %s (%s)", node.tax_id, node.name)
                    if log_taxid:
                        log_taxid(node.tax_id, node.name, 0, 0, 0)
                    continue
                elif len(seqs) < a.min_seqs_for_filtering:
                    logging.debug('%d sequence(s) for %s (%s) [action: %s]',
                                 len(seqs), node.tax_id, node.name,
                                 a.rare_taxon_action)
                    if a.rare_taxon_action == DROP:
                        continue
                    elif a.rare_taxon_action == KEEP:
                        kept_ids |= seqs
                    else:
                        raise ValueError("Unknown action: {0}".format(
                            a.rare_taxon_action))

                f = executor.submit(filter_worker,
                        sequence_file=a.sequence_file,
                        node=node,
                        seqs=seqs,
                        distance_cutoff=a.distance_cutoff,
                        log_taxid=log_taxid)
                futs[f] = {'n_seqs': len(seqs), 'node': node}

            complete = 0
            while futs:
                done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
                complete += len(done)
                sys.stderr.write('{0:8d}/{1:8d} taxa completed\r'.format(complete,
                    complete + len(pending)))
                for f in done:
                    if f.exception():
                        logging.exception("Error in child process: %s", f.exception())
                        executor.shutdown(False)
                        raise f.exception()

                    info = futs.pop(f)
                    kept = f.result()
                    kept_ids |= kept
                    if len(kept) == 0:
                        logging.info('Pruned all %d sequences for %s (%s)',
                                info['n_seqs'], info['node'].tax_id,
                                info['node'].name)
                    elif len(kept) != info['n_seqs']:
                        logging.info('Pruned %d/%d sequences for %s (%s)',
                                info['n_seqs'] - len(kept), info['n_seqs'],
                                info['node'].tax_id, info['node'].name)

        # Extract all of the sequences that passed.
        logging.info('Extracting %d sequences', len(kept_ids))
        wrap.esl_sfetch(a.sequence_file, kept_ids, fp)

    # Filter seqinfo to sequences that passed.
    if a.filtered_seqinfo:
        with open(a.seqinfo_file.name) as fp:
            r = csv.DictReader(fp)
            rows = (i for i in r if i['seqname'] in kept_ids)
            with a.filtered_seqinfo as ofp:
                w = csv.DictWriter(ofp, r.fieldnames, quoting=csv.QUOTE_NONNUMERIC)
                w.writeheader()
                w.writerows(rows)
