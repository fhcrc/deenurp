commands = ('add_reps',
            'hrefpkg_build',
            'search_sequences',
            'select_references',
            'fill_lonely',
            'tax2tree',
            'filter_outliers',
            'rdp_extract_genbank',
            'rdp_sequence_filter',
            'cluster_refs',
            'expand_named')

def itermodules(root=__name__):
    for command in sorted(commands):
        yield command.replace('_', '-'), \
                __import__('.'.join((root, command)), fromlist=[command])
