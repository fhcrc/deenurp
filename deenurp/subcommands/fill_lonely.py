"""
Find selected taxa which are the sole descendent of their parent, add
additional 'nearby' taxa.
"""

import argparse
import copy
import csv
import logging
import multiprocessing
import shutil
import sys

from Bio import SeqIO
from concurrent import futures
from taxtastic import taxtable

from .. import util, wrap

RANK = 'species'
PARENT_RANK = 'genus'

def is_lonely(node, parent_rank=PARENT_RANK):
    try:
        parent = node.at_rank(parent_rank)
        return sum(True for i in parent if i.rank == node.rank) == 1
    except KeyError:
        return False

def fill_lonely(node_id, parent_id, full_taxonomy, full_fasta):
    parent_node = full_taxonomy.get_node(parent_id)
    lonely_node = next(i for i in parent_node if i.tax_id == node_id)
    other_nodes = [i for i in parent_node if i.rank == lonely_node.rank and i != lonely_node]
    other_sequence_ids = [s for i in other_nodes for s in i.subtree_sequence_ids()]
    if other_sequence_ids:
        with util.ntf(suffix='.fasta') as tf, util.ntf(suffix='.fast.tre') as tree_fp:
            wrap.esl_sfetch(full_fasta, other_sequence_ids, tf)
            tf.seek(0)
            sequences = SeqIO.parse(tf, 'fasta')

            # Align
            sys.stderr.write('Node {0}: cmalign {1} sequences\r'.format(node_id, len(other_sequence_ids)))
            aligned = list(wrap.cmalign(sequences))

            # Run FastTree
            sys.stderr.write('Node {0}: FastTree {1} sequences\r'.format(node_id, len(other_sequence_ids)))
            wrap.fasttree(aligned, '/dev/null', tree_fp, gtr=True)
            tree_fp.close()

            # Select reps
            sys.stderr.write('Node {0}: Minimizing ADCL\r'.format(node_id))
            prune_leaves = wrap.rppr_min_adcl_tree(tree_fp.name, 5)
            return frozenset(other_sequence_ids) - frozenset(prune_leaves)

    return frozenset()

def build_parser(p):
    p.add_argument('search_fasta', help="""Sequence corpus""")
    p.add_argument('search_seqinfo', help="""Sequence corpus metadata""", type=argparse.FileType('r'))
    p.add_argument('search_taxtable', help="""Metadata for sequence corpus""", type=argparse.FileType('r'))

    p.add_argument('chosen_fasta', help="""Chosen reference sequences to augment""")
    p.add_argument('chosen_seqinfo', help="""Chosen reference sequence
            metadata""", type=argparse.FileType('r'))

    p.add_argument('output', help="Output file (fasta)", type=argparse.FileType('w'))
    p.add_argument('output_seqinfo', help="""Destination to write seqinfo for
            new representatives""", type=argparse.FileType('w'))

    rank_group = p.add_argument_group('Taxonomic Ranks')
    rank_group.add_argument('--lonely-rank', default='species', help="""Rank to
            search for lonely taxa [default: %(default)s]""")
    rank_group.add_argument('--parent-rank', default='genus', help="""Rank
            above `--lonely-rank` [default: %(default)s]""")

    thread_group = p.add_argument_group('Threading')
    thread_group.add_argument('--threads', default=multiprocessing.cpu_count(), type=int,
            help="""Number of threads [default: %(default)s]""")

def action(args):
    logging.info("Loading taxtable")
    with args.search_taxtable as fp:
        full_taxonomy = taxtable.read(fp)

    logging.info("Loading chosen sequence metadata")
    chosen_taxonomy = copy.deepcopy(full_taxonomy)
    chosen_taxonomy.populate_from_seqinfo(args.chosen_seqinfo)
    chosen_taxonomy.prune_unrepresented()

    logging.info("loading full sequence metadata")
    full_taxonomy.populate_from_seqinfo(args.search_seqinfo)

    # Find lonely
    nodes = [i for i in chosen_taxonomy if i.rank == args.lonely_rank]
    lonely_nodes = [i for i in nodes if is_lonely(i)]
    additional_reps = set()
    futs = []
    with futures.ThreadPoolExecutor(args.threads) as executor:
        for node in lonely_nodes:
            futs.append(executor.submit(fill_lonely, node.tax_id, node.at_rank(args.parent_rank).tax_id, full_taxonomy, args.search_fasta))
        while futs:
            try:
                done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
                futs = set(pending)
                for f in done:
                    if f.exception():
                        raise f.exception()
                    additional_reps |= f.result()
                sys.stderr.write("{0:6d}/{1:6d} complete        \r".format(len(lonely_nodes) - len(pending),
                        len(lonely_nodes)))
            except futures.TimeoutError:
                pass # Keep waiting
            except:
                logging.exception("Caught error in child thread - exiting")
                executor.shutdown(False)
                raise

    logging.info("%d additional references", len(additional_reps))
    with open(args.chosen_fasta) as fp, args.output as ofp:
        shutil.copyfileobj(fp, ofp)
        wrap.esl_sfetch(args.search_fasta, additional_reps, ofp)

    with args.chosen_seqinfo as fp, args.output_seqinfo as ofp, \
            args.search_seqinfo as sub_fp:
        fp.seek(0)
        r = csv.DictReader(fp)
        w = csv.DictWriter(ofp, r.fieldnames, quoting=csv.QUOTE_NONNUMERIC, lineterminator='\n')
        w.writeheader()
        w.writerows(r)

        args.search_seqinfo.seek(0)
        for row in csv.DictReader(sub_fp):
            if row['seqname'] in additional_reps:
                w.writerow(row)
