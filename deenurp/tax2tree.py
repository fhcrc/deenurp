import csv
import logging
import os

from t2t import nlevel

from . import wrap
from .tax import TaxNode


TAX2TREE_RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                  'species']

def generate_tax2tree_map(refpkg, output_fp):
    """
    Generate a tax2tree map from a reference package, writing to output_fp
    """
    with open(refpkg.file_abspath('taxonomy')) as fp:
        tax_root = TaxNode.from_taxtable(fp)

    def lineage(tax_id):
        l = tax_root.get_node(tax_id).lineage()
        d = {i.rank: i.tax_id for i in l}
        r = ('{0}__{1}'.format(rank[0], d.get(rank, '')) for rank in TAX2TREE_RANKS)
        return '; '.join(r)

    with open(refpkg.file_abspath('seq_info')) as fp:
        reader = csv.DictReader(fp)
        seq_map = ((i['seqname'], i['tax_id']) for i in reader)

        for seqname, tax_id in seq_map:
            if not tax_id:
                l = 'Unclassified'
            else:
                l = lineage(tax_id)

            print >> output_fp, '\t'.join((seqname, l))

def parse_tax2tree_out(fp):
    """
    Parse the output of tax2tree, returning the most specific taxon, or None
    """
    for line in fp:
        sequence, lineage = line.strip().split('\t')
        lineage = [i.split('__')[1] or None for i in lineage.split('; ')]
        lineage = filter(None, lineage)
        yield sequence, lineage[-1] if lineage else None

def update_taxids(refpkg, tax2tree_dict, output_fp):
    with open(refpkg.file_abspath('taxonomy')) as fp:
        tax_root = TaxNode.from_taxtable(fp)
    def lineage_ids(tax_id):
        if not tax_id:
            return frozenset()
        n = tax_root.get_node(tax_id)
        s = frozenset(i.tax_id for i in n.lineage())
        return s

    with open(refpkg.file_abspath('seq_info')) as fp:
        dialect = csv.Sniffer().sniff(fp.read(500))
        fp.seek(0)
        r = csv.DictReader(fp, dialect=dialect)
        h = r.fieldnames
        h.append('orig_tax_id')
        w = csv.DictWriter(output_fp, h, dialect=dialect)
        w.writeheader()
        for i in r:
            i['orig_tax_id'] = i['tax_id']
            n = tax2tree_dict.get(i['seqname'])
            if n and n not in lineage_ids(i['tax_id']):
                node = tax_root.get_node(n)
                new_name, new_rank = node.name, node.rank
                if i['tax_id']:
                    node = tax_root.get_node(i['tax_id'])
                    orig_name, orig_rank = node.name, node.rank
                else:
                    orig_name, orig_rank = '', ''

                logging.info('%s changed from "%s" (%s) to "%s" (%s)',
                        i['seqname'], orig_name, orig_rank, new_name, new_rank)

                i['tax_id'] = n
            w.writerow(i)

def tax2tree(refpkg, output_fp):
    """
    Run tax2tree on a reference package lacking names
    """
    with wrap.ntf() as tf, wrap.ntf() as output:
        output.close()
        consensus = output.name + '-consensus-strings'
        generate_tax2tree_map(refpkg, tf)
        tf.close()
        nlevel.main(['nlevel', '--consensus-map', tf.name,
            '--output', output.name,
            '--tree', refpkg.file_abspath('tree')])
        try:
            with open(consensus) as fp:
                d = dict(parse_tax2tree_out(fp))
                update_taxids(refpkg, d, output_fp)
        finally:
            os.remove(consensus)
