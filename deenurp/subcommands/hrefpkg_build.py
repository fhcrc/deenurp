"""
Build a hierarchical set of reference packages
"""

import csv
import random

from Bio import SeqIO
from taxtastic.refpkg import Refpkg

from .. import wrap, tax

def fasta_db(s):
    c = s.count(':')
    if not c:
        return s, None
    elif c == 1:
        return s.split(':')
    else:
        raise ValueError("invalid specification: {0}".format(s))

def build_parser(p):
    p.add_argument('sequence_file', help="""All sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""")
    p.add_argument('taxonomy', help="""Taxtable""")
    p.add_argument('--index-rank', help="""Rank for individual reference
            packages [default: %(default)s]""", default='order')
    p.add_argument('--threads', type=int, default=12, help="""Number of threads
            [default: %(default)d]""")

def action(a):
    index_rp = build_index_refpkg(a.sequence_file, a.seqinfo_file, a.taxonomy)

    with open(index_rp.file_abspath('taxonomy')) as fp:
        taxonomy = tax.TaxNode.from_taxtable(fp)

    nodes = [i for i in taxonomy if i.rank == a.index_rank]
    for i, node in enumerate(nodes):
        print '{0}: {1} ({2}/{3})'.format(node.tax_id, node.name, i + 1,
                len(nodes))
        tax_id_refpkg(index_rp, node.tax_id)

def seqinfo_tax_ids(seqinfo_fp):
    """
    Return all unique tax_ids in seqinfo_fp
    """
    r = csv.DictReader(seqinfo_fp)
    return frozenset(i['tax_id'] for i in r)

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

def choose_sequence_ids(taxonomy, seqinfo_rows, per_species=5):
    """
    Select sequences
    """
    for i in seqinfo_rows:
        taxonomy.get_node(i['tax_id']).sequence_ids.append(i['seqname'])

    species = (i for i in taxonomy if i.rank == 'species')
    for s in species:
        spec_seqs = list(s.subtree_sequence_ids())
        if len(spec_seqs) > per_species:
            spec_seqs = random.sample(spec_seqs, per_species)
        for i in spec_seqs:
            yield i

def tax_id_refpkg(index_refpkg, tax_id, threads=12):
    """
    Build a reference package containing all descendants of tax_id from an
    index reference package.
    """
    with wrap.ntf(prefix='taxonomy', suffix='.csv') as tax_fp, \
         wrap.ntf(prefix='aln_sto', suffix='.sto') as sto_fp, \
         wrap.ntf(prefix='tree', suffix='.tre') as tree_fp, \
         wrap.ntf(prefix='tree', suffix='.stats') as stats_fp, \
         wrap.ntf(prefix='seq_info', suffix='.csv') as seq_info_fp:
        stats_fp.close()

        # Subset taxonomy
        with open(index_refpkg.file_abspath('taxonomy')) as fp:
            full_tax = tax.TaxNode.from_taxtable(fp)
            descendants = set(i.tax_id for i in full_tax.get_node(tax_id))
            assert descendants
            fp.seek(0)
            r = csv.DictReader(fp)
            w = csv.DictWriter(tax_fp, r.fieldnames)
            w.writeheader()
            rows = (i for i in r if i['tax_id'] in descendants)
            w.writerows(rows)
            tax_fp.close()

        # Subset seq_info
        with open(index_refpkg.file_abspath('seq_info')) as fp:
            r = csv.DictReader(fp)
            w = csv.DictWriter(seq_info_fp, r.fieldnames, quoting=csv.QUOTE_NONNUMERIC)
            w.writeheader()
            rows = [i for i in r if i['tax_id'] in descendants]
            keep_seq_ids = frozenset(choose_sequence_ids(full_tax, rows))
            rows = [i for i in rows if i['seqname'] in keep_seq_ids]
            assert len(rows) == len(keep_seq_ids)
            w.writerows(rows)
            seq_info_fp.close()

        # Align
        sequences = SeqIO.parse(index_refpkg.file_abspath('aln_fasta'), 'fasta')
        sequences = [i for i in sequences if i.id in keep_seq_ids]
        print 'Tax id {0}: {1} sequences'.format(tax_id, len(sequences))

        # No sense in building with one sequence
        if len(sequences) == 1:
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
        rp.update_phylo_model('FastTree', stats_fp.name)
        rp.update_file('profile', wrap.CM)
        rp.commit_transaction()

        return rp
