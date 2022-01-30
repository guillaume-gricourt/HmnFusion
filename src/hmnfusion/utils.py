import json
from typing import Dict, List


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


def update_list(li: List[object], indexes:List[int]) -> List[object]:
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
