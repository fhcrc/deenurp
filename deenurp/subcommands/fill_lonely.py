"""
Add taxa of the same rank to sole descendents of their parent taxon.
"""

import argparse
import copy
import csv
import io
import logging
import shutil
import sys

from Bio import SeqIO
from concurrent import futures
from taxtastic import taxtable

from .. import config, util, wrap

RANK = 'species'
PARENT_RANK = 'genus'


def is_lonely(node, parent_rank=PARENT_RANK):
    """
    Returns whether a node is lonely [the only node with rank ``node.rank``
    under the parent].

    Nodes without a parent at ``parent_rank`` are defined to be *non-lonely*
    """
    try:
        parent = node.at_rank(parent_rank)
        return sum(True for i in parent if i.rank == node.rank) == 1
    except ValueError:
        return False


def fill_lonely_worker(
        node_id, parent_id, full_taxonomy, full_fasta, fa_idx, n_reps=5):
    """
    Finds some company for lonely taxonomic node identified by ``node_id``

    :param tax_id: Tax ID of lonely node
    :param parent_id: Tax ID of parent node (at desired rank)
    :param full_taxonomy: complete taxtable.TaxNode, populated with available
    sequences in ``full_fasta``
    :param full_fasta: Path to FASTA file with candidate references
    :returns: Sequence IDs to keep
    """
    parent_node = full_taxonomy.get_node(parent_id)
    lonely_node = next((i for i in parent_node if i.tax_id == node_id), None)
    if lonely_node is None:
        raise KeyError("Unable to find node with tax_id {0}".format(node_id))

    # Find other nodes with equivalent rank
    other_nodes = [i for i in parent_node if i.rank ==
                   lonely_node.rank and i != lonely_node]
    other_sequence_ids = [
        s for i in other_nodes for s in i.subtree_sequence_ids()]

    if len(other_sequence_ids) <= n_reps:
        return frozenset(other_sequence_ids)

    # Choose representatives
    with util.ntf(suffix='.fasta') as tf, \
            util.ntf(suffix='.fast.tre') as tree_fp:
        wrap.esl_sfetch(full_fasta, other_sequence_ids, tf, fa_idx)
        tf.seek(0)
        sequences = SeqIO.parse(io.TextIOWrapper(tf, encoding='utf8'), 'fasta')

        # Align
        sys.stderr.write(
            'Node {0}: cmalign {1} sequences\r'.format(
                node_id, len(other_sequence_ids)))
        aligned = list(wrap.cmalign(sequences))

        # Run FastTree
        sys.stderr.write('Node {0}: FastTree {1} sequences\r'.format(
            node_id, len(other_sequence_ids)))
        wrap.fasttree(aligned, tree_fp.name, gtr=True)
        tree_fp.close()

        # Select reps
        sys.stderr.write('Node {0}: Minimizing ADCL\r'.format(node_id))
        prune_leaves = wrap.rppr_min_adcl_tree(tree_fp.name, 5)
        return frozenset(other_sequence_ids) - frozenset(prune_leaves)


def build_parser(p):
    p.add_argument('search_fasta',
                   help="""Sequence corpus""")
    p.add_argument('search_seqinfo',
                   help="""Sequence corpus metadata""",
                   type=argparse.FileType('r'))
    p.add_argument('search_taxtable',
                   help="""Metadata for sequence corpus""",
                   type=argparse.FileType('r'))

    p.add_argument('chosen_fasta',
                   help="""Chosen reference sequences to augment""")
    p.add_argument('chosen_seqinfo',
                   help="""Chosen reference sequence metadata""",
                   type=argparse.FileType('r'))

    p.add_argument('output',
                   help="""Output file (fasta)""",
                   type=argparse.FileType('wb'))
    p.add_argument('output_seqinfo',
                   help="""Destination to write seqinfo
                           for new representatives""",
                   type=argparse.FileType('w'))

    rank_group = p.add_argument_group('Taxonomic Ranks')
    rank_group.add_argument('--lonely-rank',
                            metavar='TAX_RANK',
                            default='species',
                            help="""Rank to search for lonely
                                    taxa [default: %(default)s]""")
    rank_group.add_argument('--parent-rank',
                            metavar='TAX_RANK',
                            default='genus',
                            help="""Rank above `--lonely-rank`
                                    [default: %(default)s]""")

    reps_group = p.add_argument_group('Representative Choice')
    reps_group.add_argument('-n', '--number-of-reps',
                            type=int, default=5, metavar='N',
                            help="""Number of additional sequences
                                    to add for each lonely taxonomic
                                    node [default: %(default)s]""")
    reps_group.add_argument('--include-taxids',
                            metavar='FILE',
                            type=argparse.FileType('r'),
                            help=('always include given tax_ids'))
    reps_group.add_argument('--exclude-taxids',
                            metavar='FILE',
                            type=argparse.FileType('r'),
                            help=('always exclude given tax_ids'))

    thread_group = p.add_argument_group('Threading')
    thread_group.add_argument(
        '--threads',
        default=config.DEFAULT_THREADS,
        type=int,
        help="""Number of threads [default: %(default)s]""")


def action(args):
    fa_idx = wrap.read_seq_file(args.search_fasta)

    logging.info("Loading taxtable")
    with args.search_taxtable as fp:
        full_taxonomy = taxtable.read(fp)

    logging.info("Loading chosen sequence metadata")
    chosen_taxonomy = copy.deepcopy(full_taxonomy)
    chosen_taxonomy.populate_from_seqinfo(args.chosen_seqinfo)
    chosen_taxonomy.prune_unrepresented()

    logging.info("loading full sequence metadata")
    full_taxonomy.populate_from_seqinfo(args.search_seqinfo)

    if args.exclude_taxids:
        for e in args.exclude_taxids:
            e = e.strip()
            logging.info('ignoring tax_id {}'.format(e))
            full_taxonomy.get_node(e).remove_subtree()

    # Find lonely
    nodes = [i for i in chosen_taxonomy if i.rank == args.lonely_rank]
    lonely_nodes = [i for i in nodes if is_lonely(i)]
    additional_reps = set()
    futs = []
    with futures.ThreadPoolExecutor(args.threads) as executor:
        for node in lonely_nodes:
            futs.append(
                executor.submit(
                    fill_lonely_worker,
                    node.tax_id,
                    node.at_rank(args.parent_rank).tax_id,
                    full_taxonomy,
                    args.search_fasta,
                    fa_idx,
                    n_reps=args.number_of_reps))

        while futs:
            try:
                done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
                futs = set(pending)
                for f in done:
                    if f.exception():
                        raise f.exception()
                    additional_reps |= f.result()
                sys.stderr.write(
                    "{0:6d}/{1:6d} complete        \r".format(
                        len(lonely_nodes) - len(pending),
                        len(lonely_nodes)))
            except futures.TimeoutError:
                pass  # Keep waiting
            except:
                logging.exception("Caught error in child thread - exiting")
                executor.shutdown(False)
                raise

    if args.include_taxids:
        for t in args.include_taxids:
            t = t.strip()
            logging.info('including tax_id {}'.format(t))
            for s in set(full_taxonomy.get_node(t).subtree_sequence_ids()):
                logging.info('sequence {}'.format(s))
                additional_reps.add(s)

    logging.info("%d additional references", len(additional_reps))
    with open(args.chosen_fasta, 'rb') as fp, args.output as ofp:
        shutil.copyfileobj(fp, ofp)
        wrap.esl_sfetch(args.search_fasta, additional_reps, ofp, fa_idx)

    with args.chosen_seqinfo as fp, args.output_seqinfo as ofp, \
            args.search_seqinfo as sub_fp:
        fp.seek(0)
        r = csv.DictReader(fp)
        w = csv.DictWriter(
            ofp,
            r.fieldnames,
            quoting=csv.QUOTE_NONNUMERIC,
            lineterminator='\n')
        w.writeheader()
        w.writerows(r)

        args.search_seqinfo.seek(0)
        for row in csv.DictReader(sub_fp):
            if row['seqname'] in additional_reps:
                w.writerow(row)
