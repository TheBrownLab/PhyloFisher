import argparse


class CustomHelpFormatter(argparse.HelpFormatter):
    """This class can be used to make chanes in the help"""

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


class myHelpFormatter(CustomHelpFormatter, argparse.RawTextHelpFormatter):
    pass