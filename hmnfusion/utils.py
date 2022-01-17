import json


def abort(parser, msg=""):
    """Abort the program"""
    parser.error(msg)


def read_json(path):
    """Read a json file"""
    with open(path) as fid:
        return json.load(fid)


def update_list(li, indexes):
    """From a list of index, delete items of a list"""
    up = 0
    for ix in indexes:
        del li[ix-up]
        up += 1
    return li
