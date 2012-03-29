import argparse
import logging
import sys

from .. import subcommands

def main(argv=sys.argv[1:]):
    action, arguments = parse_arguments(argv)

    loglevel = {
        0: logging.WARN,
        1: logging.INFO,
        2: logging.DEBUG,
    }.get(arguments.verbosity, logging.DEBUG)

    if arguments.verbosity > 1:
        logformat = '%(levelname)s %(module)s %(lineno)s %(message)s'
    else:
        logformat = '%(message)s'

    # set up logging
    logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)

    return action(arguments)

def parse_arguments(argv):
    """
    Extract command-line arguments for different actions.
    """
    parser = argparse.ArgumentParser(description="deenurp pipeline")

    parser.add_argument('-v', '--verbose', dest='verbosity',
            action='count', default=1,
            help="Be more verbose. Specify -vv or -vvv for even more")
    parser.add_argument('-q', '--quiet', action='store_const', const=0,
            dest='verbosity', help="Suppress output")

    # Subparsers
    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_help = subparsers.add_parser('help',
            help='Detailed help for actions using help <action>')

    parser_help.add_argument('action')

    # Add actions
    actions = {}
    for name, mod in subcommands.itermodules():
        subparser = subparsers.add_parser(name, help=mod.__doc__,
                description=mod.__doc__)
        mod.build_parser(subparser)
        actions[name] = mod.action

    arguments = parser.parse_args(argv)
    action = arguments.subparser_name

    if action == 'help':
        return parse_arguments([str(arguments.action), '-h'])

    return actions[action], arguments
