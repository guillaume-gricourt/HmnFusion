import logging
import os
import shutil
from typing import List

import numpy as np
import pandas as pd
import pysam
from hmnfusion import bed as ibed
from hmnfusion import graph, region


# Helper functions
def cigar2position(cigars: List[List[int]], start: int) -> dict:
    """Construct from a cigar and a position, a position/operation.

    Parameters
    ----------
    cigars:
        A cigar coming from pysam
    start:
        A genomic coordinate

    Returns
    -------
    dict
        A dictionary with:
          - key: genomic coordinate
          - value: cigar attribute (H, S, M, ...)

    """
    data = {}
    if cigars[0][0] in [4, 5]:
        for i in range(cigars[0][1]):
            i += 1
            data[start - i] = cigars[0][0]
        del cigars[0]
    for cigar in cigars:
        op, nb = cigar[0], cigar[1]
        if op in [4, 5]:
            for i in range(nb):
                data[start + i] = op
                i += 1
        elif op == 0:
            for i in range(nb):
                data[start] = op
                start += 1
        elif op == 1:
            pass
        else:
            start += 1
    return data


def run(params: dict, bed: ibed.Bed, g: graph.Graph) -> None:
    """Main function to quantify fusion, by updating Graph object

    Parameters
    ----------
    params: dict
        Some parameters
    bed: hmnfusion.bed.Bed
        A Bed object
    g: hmnfusion.graph.Graph
        An object containing fusion data

    Return
    ------
    None
    """
    alignment = pysam.AlignmentFile(
        params["falignment"]["path"], params["falignment"]["mode"]
    )

    nodes = g.graph.nodes
    for n in nodes:
        g.graph.nodes[n]["is_skip"] = False

        if not g.graph.nodes[n]["is_interest"]:
            continue
        # Check fusion against bed.
        sub_first = pd.DataFrame(columns=ibed.Bed.HEADER)
        sub_second = pd.DataFrame(columns=ibed.Bed.HEADER)

        if g.graph.nodes[n]["fusion"].first.is_init():
            sel = bed.df.apply(
                ibed.Bed.select_bed, axis=1, args=(g.graph.nodes[n]["fusion"].first,)
            )
            sub_first = bed.df[sel]
        if g.graph.nodes[n]["fusion"].second.is_init():
            sel = bed.df.apply(
                ibed.Bed.select_bed, axis=1, args=(g.graph.nodes[n]["fusion"].second,)
            )
            sub_second = bed.df[sel]
        if len(sub_first) > 1 or len(sub_second) > 1:
            logging.warning(
                "Fusion %s is found multiple times in bed -> skipping"
                % (g.graph.nodes[n]["fusion"],)
            )
            g.graph.nodes[n]["is_skip"] = True
        if len(sub_first) + len(sub_second) == 2:
            logging.warning(
                "Fusion %s is found on left and right of breakpoint in the bed -> skipping"
                % (g.graph.nodes[n]["fusion"],)
            )
            g.graph.nodes[n]["is_skip"] = True
        if len(sub_first) + len(sub_second) == 0:
            logging.warning(
                'Fusion %s isn"t found on left or right of breakpoint in the bed -> skipping'
                % (g.graph.nodes[n]["fusion"],)
            )
            g.graph.nodes[n]["is_skip"] = True

        if g.graph.nodes[n]["is_skip"]:
            continue

        # Init.
        bed_sel = pd.DataFrame(columns=ibed.Bed.HEADER)
        r = region.Region()
        if len(sub_first) == 1:
            bed_sel = sub_first
            r = g.graph.nodes[n]["fusion"].first
        elif len(sub_second) == 1:
            bed_sel = sub_second
            r = g.graph.nodes[n]["fusion"].second
            # Swap.
            g.graph.nodes[n]["fusion"].swap_region()
        else:
            logging.warning(
                "Fusion %s, something bad happened -> skipping"
                % (g.graph.nodes[n]["fusion"],)
            )

        # Run.
        for aligned_segment in alignment.fetch(
            bed_sel.iloc[0, 0], bed_sel.iloc[0, 1], bed_sel.iloc[0, 2]
        ):
            # Filtering.
            if (
                aligned_segment.is_unmapped
                or aligned_segment.is_duplicate
                or aligned_segment.is_supplementary
            ):
                continue

            cigar2pos = cigar2position(
                aligned_segment.cigartuples, aligned_segment.reference_start
            )
            if r.position not in cigar2pos.keys():
                continue

            g.graph.nodes[n]["fusion"].evidence.depth += 1
            # Count split reads.
            if aligned_segment.has_tag("SA"):
                g.graph.nodes[n]["fusion"].evidence.split += 1
                continue

            # Count other Chrom.
            if aligned_segment.is_paired:
                if (
                    not aligned_segment.mate_is_unmapped
                    and not aligned_segment.is_unmapped
                ):
                    if (
                        aligned_segment.next_reference_id
                        != aligned_segment.reference_id
                    ):
                        g.graph.nodes[n]["fusion"].evidence.mate += 1
                        continue

            # Count reads clipped.
            count_clipped = np.zeros((2, params["clipped"]["interval"]))
            for i in range(params["clipped"]["interval"]):
                if cigar2pos.get(r.position - i - 1, 0) in [4, 5]:
                    count_clipped[0][i] = 1
                if cigar2pos.get(r.position + i + 1, 0) in [4, 5]:
                    count_clipped[1][i] = 1

            if np.max(np.sum(count_clipped, axis=1)) >= params["clipped"]["count"]:
                g.graph.nodes[n]["fusion"].evidence.clipped += 1


def write(filename: str, name: str, g: graph.Graph) -> None:
    """Write a vcf file from a list of Fusion
    Parameters
    ----------
    filename: str
        A filename to write fusion
    name: str
        Name of sample
    g: hmnfusion.graph.Graph
        A graph object to extract data

    Return
    ------
    None
    """
    # Header.
    shutil.copyfile(
        src=os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "templates",
            "vcf",
            "vcf.header.4-2.txt",
        ),
        dst=filename,
    )

    # Fusions.
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    columns.append(name)
    df = pd.DataFrame(columns=columns)

    nodes = g.graph.nodes
    for n in nodes:
        logging.debug(g.graph.nodes[n]["fusion"])
        # First.
        ident = g.graph.nodes[n]["fusion"].get_name()
        ident_1 = ident
        if g.graph.nodes[n]["fusion"].second.is_init():
            ident_1 += "-1"
        infos = ["SVTYPE=FUS"]
        infos += ["SOFT=%s" % (g.graph.nodes[n]["fusion"].get_software(),)]
        infos += ["FROM=%s" % ("-".join(g.label_build_from(n)),)]
        infos += ["CONS=%s" % (g.graph.nodes[n]["is_consensus"],)]
        infos += ["VAF=%s" % (g.graph.nodes[n]["fusion"].evidence.get_vaf(),)]
        infos += ["DP=%s" % (g.graph.nodes[n]["fusion"].evidence.depth,)]
        infos += ["SU=%s" % (g.graph.nodes[n]["fusion"].evidence.get_sum(),)]
        infos += ["SR=%s" % (g.graph.nodes[n]["fusion"].evidence.split,)]
        infos += ["PE=%s" % (g.graph.nodes[n]["fusion"].evidence.mate,)]
        infos += ["SC=%s" % (g.graph.nodes[n]["fusion"].evidence.clipped,)]

        sinfos = ":".join(infos)
        values = [
            g.graph.nodes[n]["fusion"].first.chrom,
            g.graph.nodes[n]["fusion"].first.position,
        ]
        values += [ident_1, "N", "<FUS>", ".", ".", sinfos]
        values += ["GT:VAF:DP:SU:SR:PE:SC"]

        infos_values = []
        for x in [
            "./.",
            g.graph.nodes[n]["fusion"].evidence.get_vaf(),
            g.graph.nodes[n]["fusion"].evidence.depth,
            g.graph.nodes[n]["fusion"].evidence.get_sum(),
            g.graph.nodes[n]["fusion"].evidence.split,
            g.graph.nodes[n]["fusion"].evidence.mate,
            g.graph.nodes[n]["fusion"].evidence.clipped,
        ]:
            infos_values.append(str(x))
        values.append(":".join(infos_values))

        df = pd.concat([df, pd.DataFrame([values], columns=columns)])

        ident_2 = ident
        if g.graph.nodes[n]["fusion"].second.is_init():
            ident_2 += "-2"
            # Second.
            sinfos = ";".join(["SVTYPE=FUS", "DP=.", "SU=."])
            values = [
                g.graph.nodes[n]["fusion"].second.chrom,
                g.graph.nodes[n]["fusion"].second.position,
                ident_2,
                "N",
                "<FUS>",
                ".",
                ".",
                sinfos,
                "GT:VAF:DP:SU:SR:PE:SC",
                "./.:.:.:.:.:.:.",
            ]
            df = pd.concat([df, pd.DataFrame([values], columns=columns)])
    df.to_csv(filename, mode="a", sep="\t", index=False)
