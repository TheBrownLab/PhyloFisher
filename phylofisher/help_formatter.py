import argparse
import textwrap
from datetime import date


version = '1.0.11'
today = date.today()


class CustomHelpFormatter(argparse.HelpFormatter):
    """
    This class changes the way the help output is displayed
    """

    def _format_action_invocation(self, action):
        # This removes metvar after short option
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(action.option_strings) + ' ' + args_string

    def _split_lines(self, text, width):
        # This adds 3 spaces before lines that wrap
        lines = text.splitlines()
        for i in range(0, len(lines)):
            if i >= 1:
                lines[i] = (3 * ' ') + lines[i]
        return lines


class MyHelpFormatter(CustomHelpFormatter, argparse.RawTextHelpFormatter):
    pass


def add_global_arguments(parser, optional, required, in_help, out_help, inp, pre_suf, out_dir):
    """
    Function to add argparse arguments used in all scripts

    input: argparse group class instance "optional"
           argparse group class instance "required"
           in_help = help message output for input
           out_help = help message output for output
    """

    # Required
    if inp is True:
        required.add_argument('-i', '--input', required=True, type=str, metavar='<in_dir>',
                              help=textwrap.dedent(f"""{in_help}"""))

    # Optional
    if out_dir is True:
        today_date = today.strftime("%b.%d.%Y")
        out_dir = f'{format(parser.prog).split(".")[0]}_out_{today_date}'
        optional.add_argument('-o', '--output', default=out_dir, type=str, metavar='<out_dir>',
                              help=textwrap.dedent(f"""{out_help}
                                      """))
    if pre_suf is True:
        optional.add_argument('-p', '--prefix', metavar='<prefix>', type=str, default='',
                              help=textwrap.dedent("""\
                              Prefix of input files
                              Default: NONE
                              Example: path/to/input/prefix*"""))
        optional.add_argument('-s', '--suffix', metavar='<suffix>', type=str, default='',
                              help=textwrap.dedent("""\
                              Suffix of input files.
                              Default: NONE
                              Example: path/to/input/*suffix
                              """))
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                          help=textwrap.dedent("""\
                                      Show this help message and exit.
                                      """))


def initialize_argparse(name: object, desc: object, usage: object) -> object:
    """
    This function initialized argparse
    """
    formatter = lambda prog: MyHelpFormatter(prog, max_help_position=100)

    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(prog=name,
                                     description=desc,
                                     usage=usage,
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent(f"""\
                                         additional information:
                                            Version: {version}
                                            GitHub: https://github.com/TheBrownLab/PhyloFisher
                                            Cite: doi:https://10.1371/journal.pbio.3001365
                                            """))

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    return parser, optional, required


def get_args(parser, optional, required, pre_suf=True, inp_dir=True, out_dir=True,
             in_help='Path to input directory'):

    out_help = ('Path to user-defined output directory\n'
                f'Default: ./{format(parser.prog).split(".")[0]}_out_<M.D.Y>')

    add_global_arguments(parser, optional, required,
                         in_help=in_help, out_help=out_help, inp=inp_dir,
                         pre_suf=pre_suf, out_dir=out_dir)

    parser._action_groups.append(optional)

    return parser.parse_args()
