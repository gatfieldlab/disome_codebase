#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""A modified argparse ArgumentParser with a detailed help function for
printing global help + each subparser's help together as one help output.
Aliases for subparsers are supported.

Usage:

parser = ParserWithDetailedHelp()
sub_aliases = {'A': ['longer_A', 'pretty_A'],
               'B': ['longer_B', 'pretty_B'],
               'C': ['bad_alias'],
               'D': ['better_D']}
parser.sub_aliases = sub_aliases
subparsers = parser.add_subparsers()
subparser_A = subparsers.add_parser('A', help='A is good',
                                    aliases=sub_aliases['A'])
subparser_B = subparsers.add_parser('B', help='B is better',
                                    aliases=sub_aliases['B'])
subparser_C = subparsers.add_parser('C', help='C is not so good',
                                    aliases=sub_aliases['C'])
subparser_D = subparsers.add_parser('D', help='D could be useful',
                                    aliases=sub_aliases['D'])
parser.add_argument('-d', '--detailed-help', action='detailed_help',
                    help='a detailed help with all subs',
                    selected_subs=['A', 'B', 'D'])

parser.parse_args(['-d'])
parser.print_detailed_help(['A', 'B', 'D'])

Output:
usage: ipython [-h] [-d]
               {A,longer_A,pretty_A,B,longer_B,pretty_B,C,bad_alias,D,better_D}
               ...

positional arguments:
  {A,longer_A,pretty_A,B,longer_B,pretty_B,C,bad_alias,D,better_D}
    A (longer_A, pretty_A)
                        A is good
    B (longer_B, pretty_B)
                        B is better
    C (bad_alias)
                        C is not so good
    D (better_D)        D could be useful

optional arguments:
  -h, --help            show this help message and exit
  -d, --detailed-help   a detailed help with all subs


Input source 'A (longer_A, pretty_A)':

usage: ipython A [-h]

optional arguments:
  -h, --help  show this help message and exit


Input source 'B (longer_B, pretty_B)':

usage: ipython B [-h]

optional arguments:
  -h, --help  show this help message and exit


Input source 'D (better_D)':

usage: ipython D [-h]

optional arguments:
  -h, --help  show this help message and exit

"""

import argparse


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


class _DetailedHelpAction(argparse._HelpAction):
    def __init__(self, option_strings, dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS, help=None, selected_subs=None):
        super().__init__(option_strings, dest, default, help)
        if selected_subs:
            self.selected_subs = selected_subs
        else:
            self.selected_subs = {}

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_detailed_help(self.selected_subs)
        parser.exit()


class ParserWithDetailedHelp(argparse.ArgumentParser):
    """A subclass of argparse.ArgumentParser to provide new print_help function
    with more detail than the argparse's built-in print_help function
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.register('action', 'detailed_help', _DetailedHelpAction)
        self.sub_aliases = {}

    def print_detailed_help(self, selected_subs):
        """For arguments parsers with SubParsersAction, this detailed help
        function, prints each subparser's help separately. Aliases are handled
        through helper class _DetailedHelpAction's extra selected_subs
        argument. A typical call to add_argument for this type of action is:

         parser_instance.add_argument(
            '-d', '--detailed-help', action='detailed_help', help='prints a '
            'detailed help screen and exits',
            selected_subs=['sub1', 'sub2'])

        """
        print(self.format_help())
        # retrieve subparsers from parser
        subparsers_actions = [
            action for action in self._actions
            if isinstance(action, argparse._SubParsersAction)]
        subparsers_action = subparsers_actions[0].choices
        for choice in selected_subs:
            subparser = subparsers_action[choice]
            if self.sub_aliases[choice]:
                aliases = ' ({})'.format(
                    ', '.join(self.sub_aliases[choice]))
            else:
                aliases = ''
            print("\nInput source '{}{}':\n".format(choice, aliases))
            print(subparser.format_help())
