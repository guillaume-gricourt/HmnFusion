import hashlib
from typing import List

import pandas as pd
import pysam
from natsort import natsorted
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.dimensions import ColumnDimension, DimensionHolder


def extract(finputs: List[str]) -> pd.DataFrame:
    """Extract event of interest from vcf input files

    Parameters
    ----------
    finputs: List[str]
        A list of path of VCF files

    Return
    ------
    pd.DataFrame
        Return a dataframe with:
            - index: region_id formatted as "<contig> <genomic coordinate> <base reference> <base alternative"
            - columns:
                - "event": str (deletion)
                - "contig": str, contig name
                - "start": int, genomic coordinate
                - "deletion": str, sequence
    """
    df = pd.DataFrame()
    for finput in finputs:
        vcf_in = pysam.VariantFile(finput)
        for record in vcf_in.fetch():
            # Check if variant is a deletion.
            if len(record.ref) > len(record.alts[0]):
                region_name = "%s %s %s %s" % (
                    record.contig,
                    record.pos,
                    record.ref,
                    record.alts[0],
                )
                region_id = hashlib.md5(region_name.encode("utf8")).hexdigest()
                # Check if variant is already seen.
                if region_id not in df.index:
                    df.at[region_id, "event"] = "deletion"
                    df.at[region_id, "contig"] = record.contig
                    df.at[region_id, "start"] = record.pos + 1
                    df.at[region_id, "deletion"] = record.ref[1:].upper()

                for sample in record.samples.keys():
                    df.at[region_id, sample] = True

    df["start"] = df["start"].astype(int)
    return df


def signatures(freference: str, df: pd.DataFrame) -> pd.DataFrame:
    """Identify MH motif from event

    Parameters
    ----------
    freference: str
        Path of the reference file (fasta format expected)
    df: pd.DataFrame
        A dataframe built from VCF files

    Return
    ------
    pd.DataFrame
        The df dataframe with a supplementary column "mmej_sequence" (str)

    See also
    --------
    extract()
    """

    def _signatures(x, freference):
        len_deletion = len(x["deletion"])
        start = int(x["start"])
        region = "%s:%s-%s" % (x["contig"], start, start + (2 * len_deletion))
        faidx = pysam.faidx(freference, region, split_lines=True)
        seq = "".join(faidx[1:]).upper()

        left = seq[:len_deletion]
        right = seq[len_deletion:]

        assert x["deletion"] == left

        mh_len, mh_seq = 0, ""
        for i in range(len_deletion + 1):
            motif = left[:i]
            if right.startswith(motif):
                mh_len = i
                mh_seq = motif
            if i > mh_len:
                break

        return mh_seq

    df["mmej_sequence"] = df.apply(_signatures, axis=1, args=(freference,))
    return df


def conclude(df: pd.DataFrame) -> pd.DataFrame:
    """Conclude about the presens of MMEJ signature

    Parameters
    ----------
    df: pd.DataFrame
        A dataframe with "mmej_sequence" column

    Return
    ------
    pd.DataFrame
        The df dataframe with a supplementary column "mmej_conclusion" (str)

    See also
    --------
    extract()
    signatures()
    """

    def _conclude(x):
        res = ""
        len_deletion = len(x["deletion"])
        len_mh = len(x["mmej_sequence"])
        if len_deletion <= 1 or len_mh <= 1:
            res = "."
        elif len_deletion == len_mh:
            res = "alignment ambiguous"
        elif len_deletion > 1 and len_mh < 5:
            res = "no clear signature"
        elif len_deletion > 1 and len_mh >= 5:
            res = "mmej signature"
        return res

    df["mmej_conclusion"] = df.apply(_conclude, axis=1)
    return df


def write(filename: str, df: pd.DataFrame) -> None:
    """Write dataframe to a file.
    Parameters
    ----------
    filename: str
        Path of an output file (xlsx, excel format)
    df: pd.DataFrame
        A dataframe to write

    Return
    ------
    None

    See also
    --------
    extract()
    signatures()
    conclude()
    """
    headers = [
        "contig",
        "start",
        "deletion",
        "event",
        "mmej_sequence",
        "mmej_conclusion",
    ]
    samples = [x for x in df.columns if x not in headers]
    df = df[headers + samples]

    # Sort values.
    idx, *_ = zip(
        *natsorted(
            zip(df.index, df.contig, df.start, df.deletion, df.event),
            key=lambda x: (x[1], x[2], x[3], x[4]),
        )
    )
    df = df.loc[list(idx)]

    # Replace values.
    df = df.fillna(".")
    df.replace({True: "o"}, inplace=True)

    # Write output.
    writer = pd.ExcelWriter(filename)
    df.to_excel(writer, index=False, sheet_name="mmej")
    ws = writer.sheets["mmej"]  # pull worksheet object

    # Adjust width.
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

    writer.save()