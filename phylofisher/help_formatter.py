import argparse
import textwrap


class CustomHelpFormatter(argparse.HelpFormatter):
    """
    This class changes the way the help ouput is displayed
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


def add_global_arguments(optional, required, in_help, out_help):
    """
    Function to add argparse arguments used in all scripts

    input: argparse group class instace "optional"
           argparse group class instace "required"
           in_help = help messge output for input
           out_help = help messge output for output
    """

    # Required
    required.add_argument('-i', '--input', required=True, type=str, metavar='<in_dir>',
                          help=textwrap.dedent(f"""{in_help}"""))

    # Optional
    optional.add_argument('-o', '--output', default="output", type=str, metavar='<out_dir>',
                          help=textwrap.dedent(f"""{out_help}
                                  """))
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


def initialize_argparse(name, desc, usage, dataset, input_meta):
    """
    This function initialized argparse
    """
    formatter = lambda prog: MyHelpFormatter(prog, max_help_position=100)

    parser = argparse.ArgumentParser(prog=name,
                                     description=desc,
                                     usage=usage,
                                     formatter_class=formatter,
                                     add_help=False,
                                     epilog=textwrap.dedent("""\
                                         additional information:
                                            Version: 0.1
                                            GitHub: https://github.com/DavidZihala/PhyloFisher
                                            Cite: 
                                            """))

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    if dataset:
        optional.add_argument('-d', '--dataset_folder', metavar='<dataset>', type=str,
                              help=textwrap.dedent("""\
                                      Path to directory containing the dataset. (ONLY if different from original one 
                                      in config.ini)
                                      """))
    if input_meta:
        optional.add_argument('-im', '--input_metadata', metavar='<in_meta.tsv>', type=str,
                              help=textwrap.dedent("""\
                                      Path to input metadata file in TSV format. (ONLY if different from original one 
                                      in config.ini)
                                      """))

    return parser, optional, required


def get_args(parser, optional, required,
             in_help='Path to input directory',
             out_help=('Path to alternative output directory\n'
                       'Default: ./output')):
    add_global_arguments(optional, required, in_help, out_help)

    parser._action_groups.append(optional)

    return parser.parse_args()
