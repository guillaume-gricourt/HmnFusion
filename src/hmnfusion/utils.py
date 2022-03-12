import json
import logging
import pysam
import re
import shutil
import subprocess
from typing import Dict, List


class ExecutableNotFound(Exception):
    """Custom class to throw an error if an executable is not found."""

    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return "ExecutableNotFound, {0}".format(self.message)
        else:
            return "ExecutableNotFound"


def abort(parser, msg: str = ""):
    """Abort the program

    Parameters
    ----------
    parser:
        The parser to use
    msg: str
        The message to throw from the parser

    Return
    ------
    """
    parser.error(msg)


def read_json(path: str) -> Dict:
    """Read a json file

    Parameters
    ----------
    path: str
        The path of the file to read

    Return
    ------
    Dict
        A dictionary.
    """
    with open(path) as fid:
        return json.load(fid)


def update_list(li: List[object], indexes: List[int]) -> List[object]:
    """From a list of index, delete items of a list

    Parameters
    ----------
    li: List[object]
        A list with differents kinds of items
    indexes: List[int]
        A list of index to delete in the previous list

    Return
    ------
    List[object]
        The li object updated
    """
    up = 0
    for ix in indexes:
        del li[ix - up]
        up += 1
    return li


def cmdline(args: List[str], logger: logging.Logger = logging.getLogger()) -> int:
    """Run a command line.

    Parameters
    ----------
    args: List[str]
        A list of argument
    logger: logging.Logger (default: logging.getLogger())
        A logger object

    Return
    ------
    int
        Return code from the comand line
    """
    ret = subprocess.run(args, capture_output=True, encoding="utf8")
    if ret.stdout is not None:
        logging.info(ret.stdout)
    if ret.stderr is not None:
        logging.warning(ret.stderr)
    return ret.returncode


def check_bam_index(path: str) -> bool:
    """Build index file for a bam file if not found.

    Parameters
    ----------
    path: str
        Path of the bam file

    Return
    ------
    bool
        True if index is present or is written, False otherwise.
    """
    alignment = pysam.AlignmentFile(path)
    try:
        alignment.check_index()
    except ValueError as e:
        pysam.index(path)
    finally:
        alignment.close()
    try:
        alignment = pysam.AlignmentFile(path, require_index=True)
        return True
    except Exception:
        return False


def find_executable(executable: str, msg: str = "") -> None:
    """Find an executable in the path and raise an error if not found.

    Parameters
    ----------
    executable: str
        Name of the executable
    msg: str (default: executable name)
        Message to throw in the error
    Return
    ------
    None
    """
    if shutil.which(executable) is None:
        if msg == "":
            msg = executable
        raise ExecutableNotFound(msg)


def validate_name_sample(name: str) -> str:
    """Validate sample name to fit in VCF file (no space allowed)

    Parameters
    ----------
    name: str
        Name to check

    Raises
    ------
    ValueError
        If the characters not allowed are in the string.

    Return
    ------
    str
        Name to check
    """
    if re.search(r"\s+", name):
        raise ValueError("Space are not allowed")
    return name
