import pandas as pd
from hmnfusion import region


class Bed(object):
    """The Bed class lets to store genomic coordinates data.

    Attributes
    ----------
    HEADER: List[str]
        Contains the header of the df attribute
    df: pd.DataFrame
        A dataframe containing three columns: chrom, start, end

    Methods
    -------
    __init__()
        Construct an empty object.
    """
    HEADER = ["chrom", "start", "end"]

    def __init__(self) -> None:
        self.df = pd.DataFrame()

    @classmethod
    def from_bed(cls, path: str) -> "Bed":
        """Build object from a bed file.

        Parameters
        ----------
        path: str
            A path.

        Return
        ------
        Bed
            A novel object.
        """
        bed = Bed()
        isHeader = 0
        with open(path) as fid:
            if "track" in fid.readline():
                isHeader = 1
        df = pd.read_csv(
            path,
            sep="\t",
            usecols=[0, 1, 2],
            names=Bed.HEADER,
            dtype={"chrom": str, "start": int, "end": int},
            skiprows=isHeader,
        )
        bed.df = df
        return bed

    @classmethod
    def select_bed(cls, x: pd.Series, r: region.Region) -> bool:
        """Select from a bed file, regions cross over an other region

        Parameters
        ----------
        x: pd.Series
            Columns coming from df attribute from a Bed object
        r: region.Region
            A region object

        Returns
        -------
        bool:
            True if the region is in the bed file, False otherwise
        """
        if x["chrom"] == r.chrom:
            if x["start"] <= r.position and x["end"] >= r.position:
                return True
        return False
