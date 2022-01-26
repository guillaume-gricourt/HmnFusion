import json
import re

import bs4
from hmnfusion import _version, fusion, graph, region, utils
from pysam import VariantFile

RE_GENEFUSE_LABEL = re.compile(r'[+-](\w+):(\d+)')
RE_GENEFUSE_TOTAL = re.compile(r'total:\s{1}(\d+),')
RE_LUMPY_CHROM = re.compile(r'\D*(\d+)')
RE_LUMPY_ALT = re.compile(r'(\w+):(\d+)')


# Genefuse.
def parse_genefuse_label(label):
    """Parse Genefuse item. Return Fusion"""
    f = fusion.Fusion("genefuse")
    # Fusion breakpoint.
    sleft, sright = label.split("___")
    m = re.search(RE_GENEFUSE_LABEL, sleft)
    if m:
        r = region.Region(
            m.group(1),
            m.group(2),
            "left"
        )
        f.set_region(r)
    m = re.search(RE_GENEFUSE_LABEL, sright)
    if m:
        r = region.Region(
            m.group(1),
            m.group(2),
            "right"
        )
        f.set_region(r)
    # Evidence.
    mev = re.search(RE_GENEFUSE_TOTAL, label)
    if mev:
        evidence = mev.group(1)
    f.evidence.raw = evidence
    return f


def read_genefuse_json(graph, filename):
    """Read Genefuse json file. Return list of Fusion"""
    data = utils.read_json(filename)
    for label in data.get("fusions", []).keys():
        if not label.lower().startswith("fusion"):
            continue
        graph.add_node(parse_genefuse_label(label))


def read_genefuse_html(graph, filename):
    """Read Genefuse html file. Return list of Fusion"""
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


def read_genefuse(graph, filename, form="json"):
    """Choose good function to read genefuse file"""
    if form == "json":
        read_genefuse_json(graph, filename)
    elif form == "html":
        read_genefuse_html(graph, filename)


# Lumpy.
def read_lumpy_vcf(graph, flumpy):
    """Read Lumpy vcf file. Return list of Fusion"""
    treats = set()
    vcf_in = VariantFile(flumpy)
    for record in vcf_in.fetch():
        # Check if variant is already seen.
        if (
            "SVTYPE" in record.info.keys() and
            record.info.get("SVTYPE") != "BND"
        ):
            continue
        ident_number, ident_paired = record.id.split("_")
        if ident_number in treats:
            continue
        treats.add(ident_number)

        # Build fusion.
        f = fusion.Fusion("lumpy")
        r = region.Region(
            record.chrom,
            int(record.pos)-1
        )
        f.set_region(r)

        alt = record.alts[0]
        m = re.search(RE_LUMPY_ALT, alt)
        alt_chrom, alt_pos = 0, 0
        if m:
            alt_chrom = m.group(1)
            alt_pos = int(m.group(2))-1
        r = region.Region(
            alt_chrom,
            alt_pos
        )
        f.set_region(r)

        evidence = 0
        if "SU" in record.info.keys():
            evidence = record.info.get("SU")[0]
        f.evidence.raw = evidence

        graph.add_node(f)


def read_lumpy(graph, filename, form="vcf"):
    """Choose good function to read genefuse file"""
    if form == "vcf":
        read_lumpy_vcf(graph, filename)


# Input/Output.
def read_hmnfusion_json(filename):
    """Read json file construct with exctractfusion command.
    Return a list of Fusion
    """
    data = utils.read_json(filename)
    return graph.Graph.from_dict(data)


def write_hmnfusion_json(
    filename,
    finputs,
    graph
):
    """Write list of fusion to a json file"""
    data = {}
    data["inputs"] = finputs
    data["software"] = dict(
        name=_version.__app_name__,
        version=_version.__version__
    )
    graph.update_graph_metadata(data)
    data = graph.to_dict()

    with open(filename, "w") as fod:
        json.dump(data, fod, cls=fusion.FusionEncoder)
