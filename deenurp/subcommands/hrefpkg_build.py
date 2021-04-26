"""
Build a hierarchical set of reference packages.

Chooses 5 sequences for each species; or next rank up if species-level
sequences are not available.
"""

import argparse
import copy
import csv
import functools
import itertools
import logging
import os.path
import random
import sys
import tempfile

from concurrent import futures

from Bio import SeqIO
from taxtastic.refpkg import Refpkg
from taxtastic.taxtable import TaxNode

from .. import config, wrap, util

PER_TAXON = 5

def comma_set(s):
    s = s.split(',')
    s = frozenset(i.strip() for i in s)
    return s

def build_parser(p):
    p.add_argument('sequence_file', help="""All sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""")
    p.add_argument('taxonomy', help="""Taxtable""")
    p.add_argument('--index-rank', help="""Rank for individual reference
            packages [default: %(default)s]""", default='order')
    p.add_argument('--threads', type=int, default=config.DEFAULT_THREADS,
                   help="""Number of threads [default: %(default)d]""")
    p.add_argument('--only', help="""List of taxids to keep""", type=comma_set)
    p.add_argument('--output-dir', default='.')
    partition_group = p.add_argument_group('Generate partitiond reference packages',
            """Generate reference packages containing partial coverage at
            a taxon. For validation.""")
    partition_group.add_argument('--partition-below-rank', help="""Rank above which to
            partition PARTITION_PROP taxa at [default: no pruning]""")
    partition_group.add_argument('--partition-rank', help="""Rank at which to partition
            PARTITION_PROP taxa""")
    partition_group.add_argument('--partition-log', type=argparse.FileType('w'),
            help="""Optional path to write partitiond tax_ids""")
    partition_group.add_argument('--partition-proportion', type=float, default=0.5,
        metavar='PARTITION_PROP', help="""Proportion of sequences to partition
        [default: %(default)f]""")
    p.add_argument('--seed', type=int, default=1)

def partition_hrefpkg(a, taxonomy):
    """
    Partition a taxonomy into two non-overlapping taxonomies.
    """
    with a.partition_log or util.nothing():
        p1, p2 = partition_taxonomy(taxonomy, a.partition_below_rank,
                a.partition_rank, a.partition_proportion,
                partition_log=a.partition_log)
    for i, p in enumerate((p1, p2)):
        i = str(i + 1)
        args = copy.deepcopy(a)
        args.partition_below_rank = None
        args.partition_rank = None
        assert not os.path.isdir(i)
        os.mkdir(i)
        with util.cd(i), util.ntf(prefix='taxonomy-') as tf:
            p.write_taxtable(tf)
            tf.close()
            args.taxonomy = tf.name
            logging.info("Building hrefpkg %s in %s", i, os.getcwd())
            action(args)

def action(a):
    fa_idx = wrap.read_seq_file(a.sequence_file)

    random.seed(a.seed)
    j = functools.partial(os.path.join, a.output_dir)
    if not os.path.isdir(j()):
        raise IOError('Does not exist: {0}'.format(j()))
    if os.path.exists(j('index.refpkg')):
        raise IOError('index.refpkg exists.')

    with open(a.taxonomy) as fp:
        logging.info('loading taxonomy')
        taxonomy = TaxNode.from_taxtable(fp)

    # If partitioning, partition with current args, return
    if a.partition_below_rank is not None or a.partition_rank is not None:
        if not a.partition_below_rank or not a.partition_rank:
            raise ValueError("--partition-below-rank and --partition-rank must be specified together")
        return partition_hrefpkg(a, taxonomy)

    with open(a.seqinfo_file) as fp:
        logging.info("loading seqinfo")
        seqinfo = load_seqinfo(fp)

    # Build an hrefpkg
    nodes = [i for i in taxonomy if i.rank == a.index_rank]
    hrefpkgs = []
    futs = {}
    with open(j('index.csv'), 'w') as fp, \
         open(j('train.fasta'), 'wb') as train_fp, \
         open(j('test.fasta'), 'wb') as test_fp, \
         futures.ThreadPoolExecutor(a.threads) as executor:
        def log_hrefpkg(tax_id):
            path = j(tax_id + '.refpkg')
            fp.write('{0},{0}.refpkg\n'.format(tax_id))
            hrefpkgs.append(path)

        for i, node in enumerate(nodes):
            if a.only and node.tax_id not in a.only:
                logging.info("Skipping %s", node.tax_id)
                continue

            if os.path.exists(j(node.tax_id + '.refpkg')):
                logging.warn("Refpkg exists: %s.refpkg. Skipping", node.tax_id)
                log_hrefpkg(node.tax_id)
                continue

            f = executor.submit(tax_id_refpkg, node.tax_id, taxonomy, seqinfo,
                                a.sequence_file, fa_idx, output_dir=a.output_dir,
                                test_file=test_fp, train_file=train_fp)
            futs[f] = node.tax_id, node.name

        while futs:
            done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
            for f in done:
                tax_id, name = futs.pop(f)
                r = f.result()
                if r:
                    logging.info(
                        'Finished refpkg for %s (%s) [%d remaining]', name, tax_id, len(pending))
                    log_hrefpkg(tax_id)
            assert len(futs) == len(pending)

        # Build index refpkg
        logging.info('Building index.refpkg')
        index_rp, sequence_ids = build_index_refpkg(hrefpkgs, a.sequence_file,
            seqinfo, taxonomy, fa_idx, dest=j('index.refpkg'),
            index_rank=a.index_rank)

        # Write unused seqs
        logging.info("Extracting unused sequences")
        seqs = (i for i in SeqIO.parse(a.sequence_file, 'fasta')
                if i.id not in sequence_ids)
        c = SeqIO.write(seqs, j('not_in_hrefpkgs.fasta'), 'fasta')
        logging.info("%d sequences not in hrefpkgs.", c)

def find_nodes(taxonomy, index_rank, want_rank='species'):
    """
    Find nodes to select sequences from, preferring want_rank-level nodes, but
    moving up a rank if no species-level nodes with sequences exist.
    """
    ranks = taxonomy.ranks
    rdict = dict(list(zip(ranks, list(range(len(ranks))))))
    assert index_rank in rdict
    assert want_rank in rdict

    def any_sequences_below(node):
        for i in node:
            if i.sequence_ids:
                return True
        return False

    def try_next(it):
        try:
            return next(it)
        except StopIteration:
            return None

    def inner(node):
        # Prefer nodes at want_rank with sequences
        if node.rank == want_rank and any_sequences_below(node):
            yield node
        else:
            nodes_below = itertools.chain.from_iterable(inner(i) for i in node.children)
            first = try_next(nodes_below)
            if first:
                # If child taxa have valid sequences, use them
                for i in itertools.chain([first], nodes_below):
                    yield i
            else:
                # If there are sequences here, and it's not a made up rank
                # ('below_genus', etc), and the rank is more specific than the
                # index rank, include sequences from the node.
                if (any_sequences_below(node) and 'below' not in node.rank and
                        rdict[index_rank] < rdict[node.rank]):
                    yield node
    return inner(taxonomy)

def load_seqinfo(seqinfo_fp):
    r = csv.DictReader(seqinfo_fp)
    return list(r)

def build_index_refpkg(hrefpkg_paths, sequence_file, seqinfo, taxonomy, fa_idx,
        dest='index.refpkg', **meta):
    """
    Build an index.refpkg from a set of hrefpkgs
    """

    # Clear taxonomy
    taxonomy = copy.deepcopy(taxonomy)
    for node in taxonomy:
        node.sequence_ids = set()

    def sequence_names(f):
        with open(f) as fp:
            r = csv.DictReader(fp)
            for i in r:
                yield i['seqname']

    hrefpkgs = (Refpkg(i, create=False) for i in hrefpkg_paths)
    seqinfo_files = (i.open_resource('seq_info') for i in hrefpkgs)

    # Add seqinfo
    for f in seqinfo_files:
        with f:
            taxonomy.populate_from_seqinfo(f)

    # Remove lineages without sequences
    taxonomy.prune_unrepresented()

    sequence_ids = frozenset(taxonomy.subtree_sequence_ids())

    with util.ntf(prefix='aln_fasta', suffix='.fasta') as tf, \
         util.ntf('w', prefix='seq_info', suffix='.csv') as seq_info_fp, \
         util.ntf('w', prefix='taxonomy', suffix='.csv') as tax_fp:
        wrap.esl_sfetch(sequence_file, sequence_ids, tf, fa_idx)
        tf.close()

        # Seqinfo file
        r = (i for i in seqinfo if i['seqname'] in sequence_ids)
        w = csv.DictWriter(seq_info_fp, list(seqinfo[0].keys()), lineterminator='\n',
                quoting=csv.QUOTE_NONNUMERIC)
        w.writeheader()
        w.writerows(r)
        seq_info_fp.close()

        taxonomy.write_taxtable(tax_fp)
        tax_fp.close()

        rp = Refpkg(dest, create=True)
        rp.start_transaction()
        rp.update_file('aln_fasta', tf.name)
        rp.update_file('seq_info', seq_info_fp.name)
        rp.update_file('taxonomy', tax_fp.name)
        rp.update_file('profile', wrap.CM)

        for k, v in list(meta.items()):
            rp.update_metadata(k, v)

        rp.commit_transaction()

    return rp, sequence_ids

def choose_sequence_ids(taxonomy, seqinfo_rows, per_taxon=PER_TAXON, index_rank='order'):
    """
    Select sequences
    """
    for i in seqinfo_rows:
        taxonomy.get_node(i['tax_id']).sequence_ids.add(i['seqname'])

    nodes = find_nodes(taxonomy, index_rank)
    for node in nodes:
        node_seqs = list(node.subtree_sequence_ids())
        if len(node_seqs) > per_taxon:
            random.shuffle(node_seqs)
        yield node_seqs[:per_taxon], node_seqs[per_taxon:]

def tax_id_refpkg(tax_id, full_tax, seqinfo, sequence_file, fa_idx,
        output_dir='.',
        index_rank='order', train_file=None, test_file=None):
    """
    Build a reference package containing all descendants of tax_id from an
    index reference package.
    """
    with util.ntf('w', prefix='taxonomy', suffix='.csv') as tax_fp, \
         util.ntf('w+', prefix='aln_sto', suffix='.sto') as sto_fp, \
         util.ntf('w', prefix='aln_fasta', suffix='.fasta') as fasta_fp, \
         util.ntf(prefix='tree', suffix='.tre') as tree_fp, \
         util.ntf(prefix='tree', suffix='.stats') as stats_fp, \
         util.ntf('w', prefix='seq_info', suffix='.csv') as seq_info_fp:

        # Subset taxonomy
        n = full_tax.get_node(tax_id)
        descendants = set(i.tax_id for i in n)
        assert descendants
        n.write_taxtable(tax_fp)
        tax_fp.close()

        # Subset seq_info
        w = csv.DictWriter(seq_info_fp, list(seqinfo[0].keys()),
                quoting=csv.QUOTE_NONNUMERIC)
        w.writeheader()
        rows = [i for i in seqinfo if i['tax_id'] in descendants]
        sinfo = {i['seqname']: i for i in rows}

        # Choose sequences, divide into train and test sets
        chosen = choose_sequence_ids(n, rows, index_rank=index_rank)
        keep_seq_ids = set()
        train_seq_ids = set()
        test_seq_ids = set()

        for keep, rest in chosen:
            keep_seq_ids |= frozenset(keep)
            l = len(rest)
            if l >= 2 * PER_TAXON:
                train_seq_ids |= frozenset(rest[:l // 2])
                test_seq_ids |= frozenset(rest[l // 2:])

        # Picked
        rows = [sinfo[i] for i in keep_seq_ids]
        w.writerows(rows)
        seq_info_fp.close()

        # Fetch sequences
        with util.ntf() as tf:
            wrap.esl_sfetch(sequence_file, keep_seq_ids, tf, fa_idx)
            tf.close()
            # reopen in text mode and read extracted sequences
            with open(tf.name) as seqfile:
                sequences = list(SeqIO.parse(seqfile, 'fasta'))

        logging.info("Tax id %s: %d sequences", tax_id, len(sequences))

        if len(set(str(i.seq) for i in sequences)) == 1:
            logging.warn("Skipping %s: only 1 unique sequence string", tax_id)
            return None

        # No sense in building with two sequence
        if len(sequences) < 3:
            logging.warn("Skipping: %d sequences.", len(sequences))
            return None

        # Extract training & test seqs
        if train_file:
            logging.info("%d training sequences", len(train_seq_ids))
            wrap.esl_sfetch(sequence_file, train_seq_ids, train_file, fa_idx)
        if test_file:
            logging.info("%d test sequences", len(test_seq_ids))
            wrap.esl_sfetch(sequence_file, test_seq_ids, test_file, fa_idx)

        # Cmalign
        aligned = wrap.cmalign(sequences, output=sto_fp)
        aligned = list(aligned)
        assert aligned

        # Tree
        wrap.fasttree(aligned, log_path=stats_fp.name, output_fp=tree_fp.name, threads=1, gtr=True)
        tree_fp.close()
        sto_fp.close()
        SeqIO.write(aligned, fasta_fp, 'fasta')
        fasta_fp.close()

        rp = Refpkg(os.path.join(output_dir, tax_id + '.refpkg'), create=True)
        rp.start_transaction()
        rp.update_file('aln_sto', sto_fp.name)
        rp.update_file('aln_fasta', fasta_fp.name)
        rp.update_file('tree', tree_fp.name)
        rp.update_file('seq_info', seq_info_fp.name)
        rp.update_file('taxonomy', tax_fp.name)
        try:
            rp.update_phylo_model('FastTree', stats_fp.name)
        except:
            print(stats_fp.read(), file=sys.stderr)
            raise
        rp.update_file('profile', wrap.CM)
        rp.commit_transaction()

        util.require_executable('rppr')
        rp.reroot()

        return rp.path

def partition_taxonomy(taxonomy, partition_below_rank, partition_rank, partition_prop, partition_log=None):
    """
    Trim a taxonomy of ``partition_prop`` proportion of the nodes with rank
    ``partition_rank`` below each node of rank ``partition_below_rank``
    """
    # Copy taxonomies for partitioning
    p1 = copy.deepcopy(taxonomy)
    p2 = copy.deepcopy(taxonomy)
    nodes = (i for i in p1 if i.rank == partition_below_rank)

    # Log
    if partition_log:
        writer = csv.writer(partition_log, lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
        writer.writerow(('parent', 'child', 'prune'))

        def w(parent, node, partition):
            writer.writerow((parent, node, partition))
    else:
        def w(*args):
            pass

    for node in nodes:
        children = [i for i in node if i.rank == partition_rank]
        child_count = len(children)
        partition_count = int(partition_prop * child_count)
        logging.info("Pruning %d/%d from %s-%s", partition_count, child_count,
                node.tax_id, node.name)
        prune = set(random.sample(list(range(len(children))), partition_count))

        # Lists of taxa to prune from the individual partitions
        p1_prune = [n.tax_id for i, n in enumerate(children) if i in prune]
        p2_prune = [n.tax_id for i, n in enumerate(children) if i not in prune]
        all_taxa = set(i.tax_id for i in children)
        for i, (to_prune, t) in enumerate([(p1_prune, p1), (p2_prune, p2)]):
            for j in to_prune:
                p = t.get_node(j)
                p.parent.remove_child(p)
            for j in all_taxa - set(to_prune):
                w(node.tax_id, j, i)

        # Checks
        for i in p1_prune:
            assert i in p2.index
            assert i not in p1.index
        for i in p2_prune:
            assert i in p1.index
            assert i not in p2.index

    return p1, p2
