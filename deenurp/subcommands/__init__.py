import logging
log = logging

commands = (
    'add_reps',
    'cluster_refs',
    'expand_named',
    'extract_genbank',
    'fill_lonely',
    'filter_outliers',
    'hrefpkg_build',
    'pairwise_distances',
    'partition-refs',
    'rdp_extract_genbank',
    'rdp_sequence_filter',
    'search_sequences',
    'select_references',
    'tax2tree',
    'transfer_names',
)


def itermodules(root=__name__):
    for command in sorted(commands):
        try:
            mod = __import__('.'.join((root, command)), fromlist=[command])
        except ImportError, e:
            log.error(e)
        else:
            yield (command.replace('_', '-'), mod)
