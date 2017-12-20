import argparse

from .argparser.diff import *

from .tools.diff import main_diff


# ###########################################################################
# Main function
def main():

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

    # create the argparser for the "diff" command
    diff = subparsers.add_parser(
        'diff',
        help='Differential analysis on large cohort, from Kallisto quantification.'
    )
    diff.set_defaults(func=main_diff)
    get_argparser_diff(diff)

    # recover arguments
    args = argparser.parse_args()

    # execute the command
    args.func(args, argparser)
