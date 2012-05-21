"""
Add examples of poorly-represented taxa using similarity search.

Each item in the taxonomy with less than a given number of exemplar sequences
is searched against a sequence database. Any sequences which match at high
identity, and are reverse matches to the same taxon are added to the sequence
database.
"""
import argparse
import csv
import logging
import os
import shutil

from romperroom import uclust
from taxtastic.taxtable import TaxNode

from .. import wrap, util

def build_parser(p):
    p.add_argument('sequence_file', help="""Named sequences""")
    p.add_argument('seqinfo_file', help="""Sequence info file""", type=argparse.FileType('r'))
    p.add_argument('taxonomy', help="""Taxtable""", type=argparse.FileType('r'))
    p.add_argument('unnamed_file', help="""Database of unnamed sequences""")
    p.add_argument('output', help="""Output base. %(metavar)s.fasta and
            %(metavar)s.seq_info.csv will be created""", metavar='output')
    p.add_argument('--rank', help="""Rank to operate on [default:
            %(default)s]""", default='species')
    p.add_argument('-m', '--min-at-rank', default=5, type=int, help="""Number
            of sequences at RANK, below which more representatives will be
            attempted to be recruited. [default: %(default)d]""")
    p.add_argument('--pct-id', help="""Percent ID to search at [default:
            %(default)f]""", default=0.99, type=float)

def find_underrepresented(tax_root, min_at_rank=5, rank='species'):
    """
    Find nodes of rank ``rank`` with less than ``min_at_rank`` sequences at or
    below them.
    """
    def inner(node):
        if node.rank == rank:
            sequences_below = list(node.subtree_sequence_ids())
            if len(sequences_below) < min_at_rank:
                yield node, sequences_below
        else:
            for child in node.children:
                for i in inner(child):
                    yield i
    return inner(tax_root)

def uclust_search(query, db, **kwargs):
    with util.ntf(prefix='uclust') as tf:
        uclust.search(db, query, tf.name, **kwargs)
        lines = (i for i in tf if i.startswith('H'))
        for i in uclust.parse_uclust_out(lines):
            yield i

def action(a):
    with a.taxonomy as fp:
        taxonomy = TaxNode.from_taxtable(fp)
    with a.seqinfo_file as fp:
        # List of sequences
        r = csv.DictReader(fp)
        current_seqs = frozenset(i['seqname'] for i in r)
        fp.seek(0)
        taxonomy.populate_from_seqinfo(fp)

    # Find sequences from underrepresented taxids to search
    underrep = find_underrepresented(taxonomy, a.min_at_rank, a.rank)

    tax_seqs = {}
    seq_group = {}
    for n, seqs in underrep:
        tax_seqs[n.tax_id] = seqs
        seq_group.update({i:n.tax_id for i in seqs})

    with util.ntf(prefix='to_expand-', suffix='.fasta') as expand_fp, \
         util.ntf(prefix='expand_hits-', suffix='.fasta') as hits_fp:
        # Extract sequences
        c = wrap.esl_sfetch(a.sequence_file, seq_group, expand_fp)
        logging.info('fetched %d sequences', c)
        expand_fp.close()

        # Search sequences against unnamed
        r = uclust_search(expand_fp.name, a.unnamed_file, pct_id=a.pct_id,
                maxaccepts=4, search_pct_id=0.9, trunclabels=True)
        hits = list(r)
        # Map from hit to group
        hit_group = {i.target_label: seq_group[i.query_label] for i in hits}

        # Extract hits
        c = wrap.esl_sfetch(a.unnamed_file, hit_group, hits_fp)
        logging.info('%d hits', c)
        hits_fp.close()

        # Search hits back against named file
        r = uclust_search(hits_fp.name, expand_fp.name, pct_id=a.pct_id,
                maxaccepts=1, search_pct_id=0.9, trunclabels=True)

        # Sequences which hit the same group
        update_hits = dict((i.query_label, seq_group[i.target_label])
                       for i in r
                       if seq_group[i.target_label] == hit_group[i.query_label])

        overlap = frozenset(update_hits) & current_seqs
        if overlap:
            logging.warn('%d sequences already present in corpus: %s',
                    len(overlap), ', '.join(overlap))

        # Add sequences
        with open(a.output + '.fasta', 'w') as ofp:
            with open(a.sequence_file) as fp:
                shutil.copyfileobj(fp, ofp)
            try:
                wrap.esl_sfetch(hits_fp.name, frozenset(update_hits) - current_seqs, ofp)
            finally:
                os.remove(hits_fp.name + '.ssi')

        # Write a new seq_info
        with open(a.output + '.seq_info.csv', 'w') as ofp, open(a.seqinfo_file.name) as sinfo:
            r = csv.DictReader(sinfo)
            fn = list(r.fieldnames) + ['inferred_tax_id']
            w = csv.DictWriter(ofp, fn, lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC)
            w.writeheader()
            w.writerows(i for i in r if i['seqname'] not in overlap)
            if 'cluster' in fn:
                rows = ({'seqname': k, 'tax_id': v, 'inferred_tax_id': 'yes', 'cluster': v}
                        for k, v in update_hits.items())
            else:
                rows = ({'seqname': k, 'tax_id': v, 'inferred_tax_id': 'yes'}
                        for k, v in update_hits.items())
            w.writerows(rows)
