import enum
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, List, Tuple

import pysam
from openpyxl.utils import get_column_letter
from openpyxl.worksheet import worksheet
from openpyxl.worksheet.dimensions import ColumnDimension, DimensionHolder


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


class EnumNoValue(enum.Enum):
    def __repr__(self):
        return "<%s.%s>" % (self.__class__.__name__, self.name)


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


def cmdline(args: List[str], show_output: bool = True) -> subprocess.CompletedProcess:
    """Run a command line.

    Parameters
    ----------
    args: List[str]
        A list of argument
    show_output: bool (default: True)
        Output command line

    Return
    ------
    subprocess.CompletedProcess
        Return result obtained with subprocess
    """
    ret = subprocess.run(args, capture_output=True, encoding="utf8")
    if show_output and ret.stdout is not None:
        logging.info(ret.stdout)
    if show_output and ret.stderr is not None:
        logging.warning(ret.stderr)
    return ret


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
    except ValueError:
        pysam.index(path)
    finally:
        alignment.close()
    try:
        alignment = pysam.AlignmentFile(path, require_index=True)
        return True
    except Exception:
        return False


def check_fasta_index(path: str) -> bool:
    """Build index file for a fasta file if not found.

    Parameters
    ----------
    path: str
        Path of the fasta file

    Return
    ------
    bool
        True if index is present or is written, False otherwise.
    """
    if os.path.isfile(path + ".fai"):
        return True
    try:
        pysam.faidx(path)
    except pysam.SamtoolsError:
        return False
    return True


def find_executable(name: str, toraise: bool = True) -> bool:
    """Find an executable in the PATH and raise an error if not found.

    Parameters
    ----------
    name: str
        Name of the executable
    toraise: bool (default: True)
        Raise an error, otherwise return False

    Raises
    ------
    ExecutableNotFound
        If name is not in the PATH

    Return
    ------
    bool
        True if executable is found
    """
    if shutil.which(name) is None:
        if toraise:
            raise ExecutableNotFound(name)
        else:
            return False
    return True


def validate_name_sample(name: str) -> bool:
    """Validate sample name to fit in VCF file (no space allowed)

    Parameters
    ----------
    name: str
        Name to check

    Return
    ------
    bool
        True if name is valid, False otherwise
    """
    if re.search(r"\s+", name):
        return False
    return True


def bam_to_fastq(path: str, compress: int = 4, threads: int = 1) -> Tuple[str, str]:
    """Convert a bam to two fastq files.

    Parameters
    ----------
    path: str
        Path to the bam file
    compress: int (default: 4)
        Level of compress fastq file, 0 is disable
    threads: int (default: 1)
        Number of threads to use

    Return
    ------
    Tuple[str, str]
        Path of the fastq files: forward & reverse
    """
    main_args = ["--threads", str(threads)]
    # Sort bam.
    tmp_sort = tempfile.NamedTemporaryFile(suffix=".bam")
    args = main_args + ["-n", "-o", tmp_sort.name, path]
    pysam.sort(*args)
    # Label file.
    suffixes = ["R{0}", "fastq"]
    if compress > 0:
        suffixes.append("gz")
    label_suffixes = "." + ".".join(suffixes)
    tmp_fq_fwd = tempfile.NamedTemporaryFile(
        suffix=label_suffixes.format("1"), delete=False
    )
    tmp_fq_rev = tempfile.NamedTemporaryFile(
        suffix=label_suffixes.format("2"), delete=False
    )
    tmp_fq_singleton = tempfile.NamedTemporaryFile(
        suffix=label_suffixes.format("3"), delete=False
    )
    tmp_fq_trash = tempfile.NamedTemporaryFile(
        suffix=label_suffixes.format("4"), delete=False
    )
    # Convert to fastq.
    args = main_args + [
        "-c",
        str(compress),
        "-1",
        tmp_fq_fwd.name,
        "-2",
        tmp_fq_rev.name,
        "-s",
        tmp_fq_singleton.name,
        "-0",
        tmp_fq_trash.name,
        "-n",
        tmp_sort.name,
    ]
    pysam.fastq(*args)
    # Clean up.
    tmp_sort.close()
    tmp_fq_singleton.close()
    tmp_fq_trash.close()

    return tmp_fq_fwd.name, tmp_fq_rev.name


def adjust_dim_worksheet(ws: worksheet.Worksheet) -> None:
    dim_holder = DimensionHolder(worksheet=ws)
    for ix, col in enumerate(ws.columns):
        col_nb = ix + 1
        lengths = []
        for cell in col:
            if cell.value is not None:
                lengths.append(len(str(cell.value)))
        length = max(lengths) + 6
        dim_holder[get_column_letter(col_nb)] = ColumnDimension(
            ws, min=col_nb, max=col_nb, width=length
        )
    ws.column_dimensions = dim_holder
