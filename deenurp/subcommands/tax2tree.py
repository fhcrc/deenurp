"""
Run tax2tree on a reference package, updating the seq_info file
"""
from taxtastic import refpkg

from .. import tax2tree, util

def build_parser(p):
    p.add_argument('refpkg', help="""Reference package""",
            metavar='refpkg', type=refpkg.Refpkg)

def action(args):
    with util.ntf(prefix='seq_info', suffix='.csv') as tf:
        tax2tree.tax2tree(args.refpkg, tf)
        tf.close()
        args.refpkg.update_file('seq_info', tf.name)
