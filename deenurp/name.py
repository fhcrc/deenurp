import collections
import csv
import logging
import os

from t2t import nlevel

from . import wrap

class Node(object):
    """
    Node in a taxonomy
    """
    def __init__(self, rank, tax_id, parent=None, sequence_ids=None, children=None, name=None):
        self.rank = rank
        self.name = name
        self.tax_id = tax_id
        self.parent = parent
        self.sequence_ids = sequence_ids or []
        self.children = children or []
        assert tax_id != ""

        if self.is_root:
            self.index = {self.tax_id: self}

    def add_child(self, child):
        child.parent = self
        child.index = self.index
        assert child.tax_id not in self.index
        self.index[child.tax_id] = child
        self.children.append(child)

    @property
    def is_leaf(self):
        return not self.children

    @property
    def has_sequences(self):
        return not self.sequence_ids

    @property
    def is_root(self):
        return self.parent is None

    def depth_first_iter(self):
        for child in self.children:
            for i in child.depth_first_iter():
                yield i
        yield self

    def subtree_sequence_ids(self):
        for node in self:
            for s in node.sequence_ids:
                yield s

    def path(self, tax_ids):
        """Get the node at the end of the path described by tax_ids"""
        assert tax_ids[0] == self.tax_id
        if len(tax_ids) == 1:
            return self

        n = tax_ids[1]
        try:
            child = next(i for i in self.children if i.tax_id == n)
        except StopIteration:
            raise KeyError(n)

        return child.path(tax_ids[1:])

    def get_node(self, tax_id):
        return self.index[tax_id]

    def lineage(self):
        """
        Returns all nodes between this node and the root, including this one.
        """
        if not self.parent:
            return [self]
        else:
            l = self.parent.lineage()
            l.append(self)
            return l

    def __repr__(self):
        return "<Node {0}:{1} [rank={2};children={3}]>".format(self.tax_id,
                self.name, self.rank, len(self.children))

    def __iter__(self):
        return self.depth_first_iter()

    @classmethod
    def of_taxtable(cls, taxtable_fp):
        """
        Generate a node from a taxtable fp
        """
        r = csv.reader(taxtable_fp)
        headers = next(r)
        rows = (collections.OrderedDict(zip(headers, i)) for i in r)

        row = next(rows)
        root = cls(rank=row['rank'], tax_id=row['tax_id'], name=row['tax_name'])
        for row in rows:
            rank, tax_id, name = [row[i] for i in ('rank', 'tax_id', 'tax_name')]
            path = filter(None, row.values()[4:])
            parent = root.path(path[:-1])
            parent.add_child(cls(rank, tax_id, name=name))

        return root

    @classmethod
    def of_taxdb(cls, con, root=None):
        cursor = con.cursor()
        if root is None:
            cursor.execute("SELECT tax_id, rank FROM nodes WHERE tax_id = parent_id")
        else:
            cursor.execute("SELECT tax_id, rank FROM nodes WHERE tax_id = ?", [root])

        tax_id, rank = cursor.fetchone()
        root = cls(rank=rank, tax_id=tax_id)

        def add_lineage(parent):
            cursor.execute("""SELECT tax_id, rank, tax_name
                    FROM nodes INNER JOIN names USING (tax_id)
                    WHERE parent_id = :1 and tax_id <> :1
                        AND names.is_primary = 1
                    """, [parent.tax_id])
            for tax_id, rank, name in cursor:
                node = cls(rank=rank, tax_id=tax_id, name=name)
                parent.add_child(node)
            for child in parent.children:
                add_lineage(child)

        add_lineage(root)
        return root

TAX2TREE_RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus',
                  'species']


def generate_tax2tree_map(refpkg, output_fp):
    """
    Generate a tax2tree map from a reference package, writing to output_fp
    """
    with open(refpkg.file_abspath('taxonomy')) as fp:
        tax_root = Node.of_taxtable(fp)

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
        tax_root = Node.of_taxtable(fp)
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
