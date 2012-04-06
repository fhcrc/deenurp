commands = 'add_reps', 'search_sequences', 'select_references', 'tax2tree'

def itermodules(root=__name__):
    for command in commands:
        yield command.replace('_', '-'), \
                __import__('.'.join((root, command)), fromlist=[command])
