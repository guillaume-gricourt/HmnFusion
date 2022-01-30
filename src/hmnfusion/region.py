from typing import Dict, Union


class Region(object):
    """Region let to store genomic location of a one partner of a fusion.

    Attributes
    ----------
    chrom : str
        The chromosome name.
    position: int
        Genomic coordinate.
    orientation: str, either "undefined", "left" or "right"
        The orientation of the fusion.

    Methods
    -------
    __init__(chrom="", position=0, orientation="undefined")
        Construct an object. All arguments are optional.
    """

    def __init__(self, chrom:str = "", position: int = 0, orientation: str = "undefined"):
        self.chrom = chrom
        self.position = position
        self.orientation = orientation

    # Getters Setters
    @property
    def chrom(self) -> str:
        return self._chrom

    @chrom.setter
    def chrom(self, chrom: str) -> None:
        self._chrom = chrom

    @property
    def position(self) -> int:
        return self._position

    @position.setter
    def position(self, position: int) -> None:
        self._position = int(position)

    @property
    def orientation(self) -> str:
        return self._orientation

    @orientation.setter
    def orientation(self, orientation: str) -> None:
        if orientation in ["left", "right"]:
            self._orientation = orientation
        else:
            self._orientation = "undefined"

    # Others
    def is_init(self) -> bool:
        """Chect if a Region is initialized.
        Only based on chromosome information

        Return
        ------
        bool
            True if the Region information is set.
        """
        if self._chrom == "":
            return False
        return True

    # Import Export
    def to_dict(self) -> Dict[str, Union[str, int]]:
        """Export object as a dict.

        Return
        ------
        Dict[str, Union[str, int]]
            Object represented by a dictionary.
        """
        return dict(
            orientation=self._orientation, chrom=self._chrom, position=self._position
        )

    @classmethod
    def from_dict(cls, data: Dict[str, Union[str, int]]) -> "Region":
        """Build object from a dictionary.

        Parameters
        ----------
        data
            A dictionary with all data needed to build an Region.

        Return
        ------
        Region
            A novel object.
        """
        region = Region()
        region.orientation = data.get("orientation", "")
        region.chrom = data.get("chrom", "")
        region.position = data.get("position", 0)
        return region

    # Meta functions
    def __key(self):
        return (self._orientation, self._chrom, self._position)

    def __repr__(self) -> str:
        return "%s %s:%s" % (self._orientation, self._chrom, self._position)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other: object) -> None:
        if isinstance(other, Region):
            return self.__key() == other.__key()
        return NotImplemented
