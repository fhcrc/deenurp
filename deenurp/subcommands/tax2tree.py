"""
Run tax2tree on a reference package, updating the seq_info file
"""
from taxtastic import refpkg

from .. import name, wrap

def build_parser(p):
    p.add_argument('refpkg', help="""Reference package""",
            metavar='refpkg', type=refpkg.Refpkg)

def action(args):
    with wrap.ntf(prefix='seq_info', suffix='.csv') as tf:
        name.tax2tree(args.refpkg, tf)
        tf.close()
        args.refpkg.update_file('seq_info', tf.name)
