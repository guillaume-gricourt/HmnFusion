from typing import Dict, Union


class Evidence(object):
    """Evidence let to store data to quantify fusion.

    Attributes
    ----------
    raw : int
        Number of pieces undifferentiated.
    split: int
        Number of split evidence.
    mate: int
        Number of mate evidence.
    clipped: int
        Number of clipped evidence.
    depth:
        Number of all pieces.

    Methods
    -------
    __init__(raw=0)
        Construct an object with a raw argument.
    """

    def __init__(self, raw: Union[str, float, int] = 0) -> None:
        self._raw = Evidence.set_number(raw)
        self._split = Evidence.set_number()
        self._mate = Evidence.set_number()
        self._clipped = Evidence.set_number()
        self._depth = Evidence.set_number()

    # Getters Setters
    @property
    def raw(self) -> int:
        return self._raw

    @raw.setter
    def raw(self, n: Union[str, float, int]) -> None:
        self._raw = Evidence.set_number(n)

    @property
    def split(self) -> int:
        return self._split

    @split.setter
    def split(self, n: Union[str, float, int]) -> None:
        self._split = Evidence.set_number(n)

    @property
    def mate(self) -> int:
        return self._mate

    @mate.setter
    def mate(self, n: Union[str, float, int]) -> None:
        self._mate = Evidence.set_number(n)

    @property
    def clipped(self) -> int:
        return self._clipped

    @clipped.setter
    def clipped(self, n: Union[str, float, int]) -> None:
        self._clipped = Evidence.set_number(n)

    @property
    def depth(self) -> int:
        return self._depth

    @depth.setter
    def depth(self, n: Union[str, float, int]) -> None:
        self._depth = Evidence.set_number(n)

    # Others
    @classmethod
    def set_number(cls, n: Union[str, float, int] = 0) -> int:
        """Get a number unsigned.

        Parameters
        ----------
        n : Union[str, float, int]
            A number to round and unsigned.

        Return
        ------
        int
            An absolute number
        """
        return abs(int(float(n)))

    def get_sum(self) -> int:
        """Get summary of all evidences.

        Return
        ------
        int
            Sum of all evidences.
        """
        return self._split + self._mate + self._clipped

    def get_vaf(self) -> str:
        """Return variable frequency allelic.
        Calculated as:
            VAF = sum of all evidence / total of pieces * 100

        Return
        ------
        str
            Number rounded with two decimals.
        """
        vaf = 0.0
        if self._depth > 0:
            vaf = (self.get_sum() / self._depth) * 100
        return "{:.2f}".format(vaf).replace(".", ",")

    def get_max_count(self) -> int:
        """Get maximum among raw attribute and sum of all evidences.

        Return
        ------
        int
            Maximum.
        """
        return max([self._raw, self.get_sum()])

    # Import Export
    def to_dict(self) -> Dict[str, int]:
        """Export object as a dict.

        Return
        ------
        Dict[str, int]
            Object represented by a dictionary.
        """
        return dict(
            raw=self._raw,
            split=self._split,
            mate=self._mate,
            clipped=self._clipped,
            depth=self._depth,
        )

    @classmethod
    def from_dict(cls, data: Dict[str, int]) -> "Evidence":
        """Build object from a dictionary.

        Parameters
        ----------
        data
            A dictionary with all data needed to build an Evidence.

        Return
        ------
        Evidence
            A novel object.
        """
        e = Evidence()
        e.raw = Evidence.set_number(data.get("raw", 0))
        e.split = Evidence.set_number(data.get("split", 0))
        e.mate = Evidence.set_number(data.get("mate", 0))
        e.clipped = Evidence.set_number(data.get("clipped", 0))
        e.depth = Evidence.set_number(data.get("depth", 0))
        return e

    # Meta functions
    def __key(self):
        return (self._raw, self._split, self._mate, self._clipped, self._depth)

    def __repr__(self):
        return "Ev " + " ".join(str(x) for x in self.__key())

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other: object):
        if isinstance(other, Evidence):
            return self.__key() == other.__key()
        return NotImplemented

    def __lt__(self, other: "Evidence"):
        return self.get_max_count() < other.get_max_count()

    def __gt__(self, other: "Evidence"):
        return self.get_max_count() > other.get_max_count()
