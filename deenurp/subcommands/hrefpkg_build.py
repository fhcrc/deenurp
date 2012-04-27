"""
Build a hierarchical set of reference packages.

Chooses 5 sequences for each species; or next rank up if species-level sequences are not available.
"""

import csv
import itertools
import logging
import os.path
import random
import sys
import tempfile


from Bio import SeqIO
from taxtastic.refpkg import Refpkg

from .. import wrap, tax

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
    p.add_argument('--threads', type=int, default=12, help="""Number of threads
            [default: %(default)d]""")
    p.add_argument('--only', help="""List of taxids to keep""", type=comma_set)
    p.add_argument('--seed', type=int, default=1)

def action(a):
    random.seed(a.seed)
    if not os.path.exists('index.refpkg'):
        index_rp = build_index_refpkg(a.sequence_file, a.seqinfo_file, a.taxonomy,
                index_rank=a.index_rank)
    else:
        logging.warn('index.refpkg exists. using.')
        index_rp = Refpkg('index.refpkg', create=False)

    with open(index_rp.file_abspath('taxonomy')) as fp:
        logging.info('loading taxonomy')
        taxonomy = tax.TaxNode.from_taxtable(fp)

    with open(a.seqinfo_file) as fp:
        logging.info("loading seqinfo")
        seqinfo = load_seqinfo(fp)

    nodes = [i for i in taxonomy if i.rank == a.index_rank]
    with open('index.csv', 'w') as fp:
        for i, node in enumerate(nodes):
            if a.only and node.tax_id not in a.only:
                logging.info("Skipping %s", node.tax_id)
                continue

            logging.info("%s: %s (%d/%d)", node.tax_id, node.name, i+1, len(nodes))
            if os.path.exists(node.tax_id + '.refpkg'):
                logging.warn("Refpkg exists: %s.refpkg. Skipping", node.tax_id)
                fp.write('{0},{0}.refpkg\n'.format(node.tax_id))
                continue
            r = tax_id_refpkg(index_rp, node.tax_id, seqinfo)
            if r:
                fp.write('{0},{0}.refpkg\n'.format(node.tax_id))

def find_nodes(taxonomy, index_rank, want_rank='species'):
    """
    Find nodes to select sequences from, preferring want_rank-level nodes, but
    moving up a rank if no species-level nodes with sequences exist.
    """
    ranks = taxonomy.ranks
    rdict = dict(zip(ranks, xrange(len(ranks))))
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
                # ('below_genus, etc), and the rank is more specific than the
                # index rank, include sequences from the node.
                if (any_sequences_below(node) and 'below' not in node.rank and
                        rdict[index_rank] < rdict[node.rank]):
                    yield node
    return inner(taxonomy)

def load_seqinfo(seqinfo_fp):
    r = csv.DictReader(seqinfo_fp)
    return list(r)

def build_index_refpkg(sequence_file, seqinfo_file, taxonomy, dest='index.refpkg', **meta):
    rp = Refpkg(dest, create=True)
    rp.start_transaction()
    rp.update_file('aln_fasta', sequence_file)
    rp.update_file('seq_info', seqinfo_file)
    rp.update_file('taxonomy', taxonomy)

    for k, v in meta.items():
        rp.update_metadata(k, v)

    rp.commit_transaction()

    return rp

def choose_sequence_ids(taxonomy, seqinfo_rows, per_taxon=5, index_rank='order'):
    """
    Select sequences
    """
    for i in seqinfo_rows:
        taxonomy.get_node(i['tax_id']).sequence_ids.append(i['seqname'])

    nodes = find_nodes(taxonomy, index_rank)
    for node in nodes:
        node_seqs = list(node.subtree_sequence_ids())
        if len(node_seqs) > per_taxon:
            node_seqs = random.sample(node_seqs, per_taxon)
        for i in node_seqs:
            yield i

def tax_id_refpkg(index_refpkg, tax_id, seqinfo, threads=12, index_rank='order'):
    """
    Build a reference package containing all descendants of tax_id from an
    index reference package.
    """
    with wrap.ntf(prefix='taxonomy', suffix='.csv') as tax_fp, \
         wrap.ntf(prefix='aln_sto', suffix='.sto') as sto_fp, \
         wrap.ntf(prefix='tree', suffix='.tre') as tree_fp, \
         wrap.ntf(prefix='tree', suffix='.stats') as stats_fp, \
         wrap.ntf(prefix='seq_info', suffix='.csv') as seq_info_fp:

        # Subset taxonomy
        with open(index_refpkg.file_abspath('taxonomy')) as fp:
            full_tax = tax.TaxNode.from_taxtable(fp)
            n = full_tax.get_node(tax_id)
            descendants = set(i.tax_id for i in n)
            assert descendants
            n.write_taxtable(tax_fp)
            tax_fp.close()

        # Subset seq_info
        w = csv.DictWriter(seq_info_fp, seqinfo[0].keys(), quoting=csv.QUOTE_NONNUMERIC)
        w.writeheader()
        rows = [i for i in seqinfo if i['tax_id'] in descendants]
        keep_seq_ids = frozenset(choose_sequence_ids(full_tax, rows,
                                 index_rank=index_rank))
        rows = [i for i in rows if i['seqname'] in keep_seq_ids]
        assert len(rows) == len(keep_seq_ids)
        w.writerows(rows)
        seq_info_fp.close()

        # Align
        sequences = SeqIO.parse(index_refpkg.file_abspath('aln_fasta'), 'fasta')
        with tempfile.NamedTemporaryFile() as tf:
            wrap.esl_sfetch(index_refpkg.file_abspath('aln_fasta'),
                            keep_seq_ids, tf)
            # Rewind
            tf.seek(0)
            sequences = list(SeqIO.parse(tf, 'fasta'))
        logging.info("Tax id %s: %d sequences", tax_id, len(sequences))

        if len(set(str(i.seq) for i in sequences)) == 1:
            logging.warn("Skipping %s: only 1 unique sequence string", tax_id)
            return None

        # No sense in building with one sequence
        if len(sequences) < 2:
            logging.warn("Skipping: %d sequences.", len(sequences))
            return None

        # Cmalign
        aligned = wrap.cmalign(sequences, output=sto_fp.name, mpi_args=['-np', str(threads)])
        # Tree
        wrap.fasttree(aligned, stats_fp.name, tree_fp, threads=threads, gtr=True)
        tree_fp.close()
        sto_fp.close()

        rp = Refpkg(tax_id + '.refpkg')
        rp.start_transaction()
        rp.update_file('aln_sto', sto_fp.name)
        rp.update_file('tree', tree_fp.name)
        rp.update_file('seq_info', seq_info_fp.name)
        rp.update_file('taxonomy', tax_fp.name)
        try:
            rp.update_phylo_model('FastTree', stats_fp.name)
        except:
            print >> sys.stderr, stats_fp.read()
            raise
        rp.update_file('profile', wrap.CM)
        rp.commit_transaction()

        return rp
