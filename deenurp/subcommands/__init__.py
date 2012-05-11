commands = ('add_reps', 'hrefpkg_build', 'search_sequences',
            'select_references', 'tax2tree', 'filter_outliers', 'cluster_refs', 'expand_named')

def itermodules(root=__name__):
    for command in commands:
        yield command.replace('_', '-'), \
                __import__('.'.join((root, command)), fromlist=[command])
