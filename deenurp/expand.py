import collections
import csv

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
