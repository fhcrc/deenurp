commands = (
    'add_reps',
    'cluster_refs',
    'expand_named',
    'fill_lonely',
    'filter_outliers',
    'hrefpkg_build',
    'pairwise_distances',
    'rdp_extract_genbank',
    'rdp_sequence_filter',
    'search_sequences',
    'select_references',
    'tax2tree',
    'transfer_names',
)


def itermodules(root=__name__):
    for command in sorted(commands):
        yield (command.replace('_', '-'),
               __import__('.'.join((root, command)), fromlist=[command]))
