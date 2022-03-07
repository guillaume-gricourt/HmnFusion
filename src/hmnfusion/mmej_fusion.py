import json
import pysam
import re
import subprocess

import networkx as nx
from hmnfusion import extractfusion, fusion, graph, region


def load_file(path:str , fmt: str, g: graph.Graph = graph.Graph()) -> graph.Graph:
    if fmt == "genefuse_html":
        extractfusion.read_genefuse_html(g, path)
    elif fmt == "genefuse_json":
        extractfusion.read_genefuse_json(g, path)
    elif fmt == "lumpy_vcf":
        extractfusion.read_lumpy_vcf(g, path)
    elif fmt == "hmnfusion_json":
        g = extractfusion.read_hmnfusion_json(path)
    else:
        return NotImplemented
    return g


def subset_graph(g: graph.Graph) -> List[fusion.Fusion]:
    nodes = [x for x in g.nodes if not g.nodes[x]["is_consensus"]]
    nodes = [x["fusion"] for x in nodes]
    return nodes


def fetch_reference(region: region.Region, interval: int, path_reference: str):
    start = region - interval
    end = region + interval
    reg = "%s:%s-%s" % (self.region_a.chr, start, end)
    seq = pysam.faidx(path_reference, reg)
    return seq


def filter_sequence(path: str, fus: fusion.Fusion, interval: int = 300) -> str:

    bam_original = pysam.AlignmentFile(path)
    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    with pysam.AlignmentFile(tmpfile.name, "wb", header=bam_original.header) as fod:
        for region in [fus.first, fus.second]:
            sregion = "%s:%s-%s" % (region.chrom, region.position - interval, region.position + interval)
            for aligned_segment in alignment.fetch(region = sregion):
                # Filtering.
                if aligned_segment.is_unmapped or aligned_segment.is_duplicate or aligned_segment.is_supplementary:
                    continue
                # Count split reads.
                if aligned_segment.has_tag("SA"):
                    fod.write(aligned_segment)
                    continue
                # Count other Chrom.
                if aligned_segment.is_paired:
                    if not aligned_segment.mate_is_unmapped and not aligned_segment.is_unmapped:
                        if aligned_segment.next_reference_id != aligned_segment.reference_id:
                            fod.write(aligned_segment)
                            continue
                # Count reads clipped.
                for cigar in aligned_segment.cigartuples:
                    if cigar[0] in [4, 5] and cigar[1] >= 6:
                        fod.write(aligned_segment)
                        continue
    return tmpfile.name


def create_consensus(path_reference: str, path_bam: str) -> str:
    tmpfile = tempfile.NamedTemporaryFile(delete=False)
    args = [
        "gencore",
        "-i", path,
        "-o", tmpfile.name,
        "-r", path_reference,
        "-s", "2"
    ]
    subprocess.run(args)

    return tmpfile.name


def extract_consensus(path: str, fus: fusion.Fusion, interval: int = 30) -> str:
    pass
