# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package
# in LICENSE.txt.

"""
    Common functions for path manipulation and use
"""

import logging
import os
from typing import List, Optional, Union


def get_full_path(current_file: str, *args: Union[str, List[str]]) -> str:
    """ Generate the absolute path of the directory containing a file
        e.g. __file__ == os.path.join(get_full_path(__file__), filename)

        Can be given extra arguments which will be added to the result
        to generate the absolute path of a file without using both
        get_full_path and os.path.join.
        e.g. get_full_path(__file__, "data", "pfam.hmm") -> "$DIR/data/pfam.hmm"

        Arguments:
            current_file: The file from which to take the base directory from
            *args: Strings to be added to the resulting directory

        Returns:
            A string containing the fully generated path
    """
    base = os.path.dirname(os.path.abspath(current_file))
    if not args:
        return base
    extra = os.path.join(*args)  # type: ignore  # mypy can't manage the splat
    assert isinstance(extra, str)
    return os.path.join(base, extra)


def locate_file(path: str, silent: bool = False) -> Optional[str]:
    """ Checks that a given file path is valid and that read permissions exist
        for the file

        Arguments:
            path: the file path to check
            silent: if True, debug messages won't be logged

        Returns:
            The path if it was valid or None if not
    """
    if os.path.split(path)[0]:
        if os.path.isfile(path) and os.access(path, os.R_OK):
            if not silent:
                logging.debug("Found file %r", path)
            return path
    return None
