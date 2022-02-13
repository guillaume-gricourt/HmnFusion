import json
import re

import bs4
from hmnfusion import _version, fusion, graph, region, utils
from pysam import VariantFile

RE_GENEFUSE_LABEL = re.compile(r"[+-](\w+):(\d+)")
RE_GENEFUSE_TOTAL = re.compile(r"total:\s{1}(\d+),")
RE_LUMPY_CHROM = re.compile(r"\D*(\d+)")
RE_LUMPY_ALT = re.compile(r"(\w+):(\d+)")


# Genefuse.
def parse_genefuse_label(label: str) -> fusion.Fusion:
    """From a Genefuse label create a Fusion.

    Parameters
    ----------
    label: str
        A label coming from Genefuse

    Return
    ------
    fusion.Fusion
        A novel fusion object
    """
    f = fusion.Fusion("genefuse")
    # Fusion breakpoint.
    sleft, sright = label.split("___")
    m = re.search(RE_GENEFUSE_LABEL, sleft)
    if m:
        r = region.Region(m.group(1), int(m.group(2)), "left")
        f.set_region(r)
    m = re.search(RE_GENEFUSE_LABEL, sright)
    if m:
        r = region.Region(m.group(1), int(m.group(2)), "right")
        f.set_region(r)
    # Evidence.
    mev = re.search(RE_GENEFUSE_TOTAL, label)
    if mev:
        evidence = mev.group(1)
    f.evidence.raw = int(evidence)
    return f


def read_genefuse_json(graph: graph.Graph, filename: str) -> None:
    """Read Genefuse JSON file.

    Parameters
    ----------
    graph: graph.Graph
        Add new fusions to this object
    filename: str
        A path from a JSON file

    Return
    ------
    None

    See also
    --------
    read_genefuse()
    read_genefuse_html()
    """
    data = utils.read_json(filename)
    for label in data.get("fusions", []).keys():
        if not label.lower().startswith("fusion"):
            continue
        graph.add_node(parse_genefuse_label(label))


def read_genefuse_html(graph: graph.Graph, filename: str) -> None:
    """Read Genefuse HTML file.

    Parameters
    ----------
    graph: graph.Graph
        Add new fusions to this object
    filename: str
        A path from an HTML file

    Return
    ------
    None

    See also
    --------
    read_genefuse()
    read_genefuse_json()
    """
    # Read.
    bhtml = ""
    with open(filename, encoding="utf-8") as fid:
        bhtml = fid.read()
    soup = bs4.BeautifulSoup(bhtml, features="lxml")

    # Parsing.
    div_menu = soup.find("div", attrs={"id": "menu"})
    for record in div_menu.findAll("a"):
        label = record.text
        if "fusion" not in label.lower():
            continue
        graph.add_node(parse_genefuse_label(label))


def read_genefuse(graph: graph.Graph, filename: str, form: str = "json") -> None:
    """Read Genefuse HTML file.

    Parameters
    ----------
    graph: graph.Graph
        Add new fusions to this object
    filename: str
        A path from a file
    form: str
        File format. Choices from "json" or "html"
    Return
    ------
    None

    See also
    --------
    read_genefuse_json()
    read_genefuse_html()
    """
    if form == "json":
        read_genefuse_json(graph, filename)
    elif form == "html":
        read_genefuse_html(graph, filename)


# Lumpy.
def read_lumpy_vcf(graph: graph.Graph, flumpy: str) -> None:
    """Read a Lumpy VCF file.

    Parameters
    ----------
    graph: graph.Graph
        Add new fusions to this object
    flumpy: str
        A path from a VCF file

    Return
    ------
    None

    See also
    --------
    read_lumpy()
    """
    treats = set()
    vcf_in = VariantFile(flumpy)
    for record in vcf_in.fetch():
        # Check if variant is already seen.
        if "SVTYPE" in record.info.keys() and record.info.get("SVTYPE") != "BND":
            continue
        ident_number, ident_paired = record.id.split("_")
        if ident_number in treats:
            continue
        treats.add(ident_number)

        # Build fusion.
        f = fusion.Fusion("lumpy")
        r = region.Region(record.chrom, int(record.pos) - 1)
        f.set_region(r)

        alt = record.alts[0]
        m = re.search(RE_LUMPY_ALT, alt)
        alt_chrom, alt_pos = "0", 0
        if m:
            alt_chrom = m.group(1)
            alt_pos = int(m.group(2)) - 1
        r = region.Region(alt_chrom, alt_pos)
        f.set_region(r)

        evidence = 0
        if "SU" in record.info.keys():
            evidence = record.info.get("SU")[0]
        f.evidence.raw = evidence

        graph.add_node(f)


def read_lumpy(graph: graph.Graph, filename: str, form: str = "vcf") -> None:
    """Read a Lumpy file.

    Parameters
    ----------
    graph: graph.Graph
        Add new fusions to this object
    filename: str
        A path from a file
    form: str
        File format. Only "vcf" is allowed.

    Return
    ------
    None

    See also
    --------
    read_lumpy_vcf()
    """
    if form == "vcf":
        read_lumpy_vcf(graph, filename)


# Input/Output.
def read_hmnfusion_json(filename: str) -> graph.Graph:
    """Read a JSON file built with the extractfusion command.

    Parameters
    -----------
    filename: str
        A path from a JSON file.

    Return
    ------
    graph.Graph
        A novel object

    See also
    --------
    write_hmnfusion_json()
    """
    data = utils.read_json(filename)
    return graph.Graph.from_dict(data)


def write_hmnfusion_json(filename: str, finputs: dict, g: graph.Graph) -> None:
    """Write fusions to a JSON file.

    Parameters
    -----------
    filename: str
        A path from a JSON file.
    finputs: dict
        Some metadata
    graph: graph.Graph
        A graph object

    Return
    ------
    None

    See also
    --------
    read_hmnfusion_json()
    """
    data = {}
    data["inputs"] = finputs
    data["software"] = dict(name=_version.__app_name__, version=_version.__version__)
    g.update_graph_metadata(data)
    data = g.to_dict()

    with open(filename, "w") as fod:
        json.dump(data, fod, cls=fusion.FusionEncoder)
