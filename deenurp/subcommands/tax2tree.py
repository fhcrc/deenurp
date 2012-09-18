"""
Run tax2tree on a reference package, updating the seq_info file
"""
from functools import partial

from taxtastic import refpkg

from .. import tax2tree, util

def build_parser(p):
    p.add_argument('refpkg', help="""Reference package""",
            metavar='refpkg', type=partial(refpkg.Refpkg, create=False))
    g = p.add_mutually_exclusive_group()
    g.add_argument('--allow-rename', action='store_true', default=True,
            help="""Allow sequences to be renamed [default: %(default)s]""")
    g.add_argument('--no-allow-rename', action='store_false', dest='allow_rename')

def action(args):
    with util.ntf(prefix='seq_info', suffix='.csv') as tf:
        tax2tree.tax2tree(args.refpkg, tf, allow_rename=args.allow_rename)
        tf.close()
        args.refpkg.update_file('seq_info', tf.name)
