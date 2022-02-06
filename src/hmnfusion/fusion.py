import copy
from json import JSONEncoder
from typing import Any, Dict

from hmnfusion import _version
from hmnfusion import evidence as hmn_evidence
from hmnfusion import region


class Fusion:
    """Fusion represents a genomic fusion location.

    Attributes
    ----------
    first : Region
        A first partner.
    second: Region
        A second partner.
    evidence: Evidence
        Couting evidence supporting a fusion.
    number: int
        The nth fusion.
    is_consensus: bool
        True if the fusion is built from others.
    software: str
        Name of the software from the fusion is built.

    Methods
    -------
    __init__(software:str = name of the app)
        Construct an object with a software argument.
    """

    def __init__(self, software: str = _version.__app_name__) -> None:
        self._first = region.Region()
        self._second = region.Region()
        self._evidence = hmn_evidence.Evidence()
        self._number = 0
        self._is_consensus = False
        self._software = software

    # Getters Setters
    @property
    def first(self) -> region.Region:
        return self._first

    @first.setter
    def first(self, first: region.Region) -> None:
        self._first = copy.deepcopy(first)

    @property
    def second(self) -> region.Region:
        return self._second

    @second.setter
    def second(self, second: region.Region) -> None:
        self._second = copy.deepcopy(second)

    @property
    def evidence(self) -> hmn_evidence.Evidence:
        return self._evidence

    @evidence.setter
    def evidence(self, e: hmn_evidence.Evidence) -> None:
        self._evidence = copy.deepcopy(e)

    @property
    def number(self) -> int:
        return self._number

    @number.setter
    def number(self, n: int) -> None:
        self._number = n

    @property
    def is_consensus(self) -> bool:
        return self._is_consensus

    @is_consensus.setter
    def is_consensus(self, b: bool) -> None:
        self._is_consensus = b

    @property
    def software(self) -> str:
        return self._software

    @software.setter
    def software(self, software: str) -> None:
        self._software = software

    # Others
    def update(self, other: "Fusion") -> None:
        """Replace first, second, evidence and software attribute by an other fusion.

        Parameters
        ----------
        other : fusion.Fusion
            An other fusion object.

        Return
        ------
        None
        """
        self.first = other.first
        self.second = other.second
        self.evidence = other.evidence
        self.software = other.software

    def get_name(self) -> str:
        """Get a formatted name of the attribute software.
        Build with:
            Three first char in uppercase and "_" then the number of the fusion.

        Return
        ------
        str
            A formatted name
        """
        n = self.software[:3].upper()
        if self.is_consensus:
            n = _version.__app_name__[:3].upper()
        return "%s_%s" % (n, self._number)

    def get_software(self) -> str:
        """Get the name of the software used to construct the fusion.

        Return
        ------
        str
            The entire name of the app.
        """
        if self.is_consensus:
            return _version.__app_name__
        return self.software

    def is_near(self, other: "Fusion", consensus_interval: int) -> bool:
        """Define if two fusions are close enough.

        Parameters
        ----------
        other: Fusion
            An another fusion
        consensus_interval: int
            A distance to evaluate if chromosomic location are close enough.
            Only one partner is needed to be close enough.

        Return
        ------
        bool
            True if the two fusions are near.
        """
        gaps = []
        if self.first.chrom == other.first.chrom:
            gaps.append(abs(self.first.position - other.first.position))
            if self.second.chrom == other.second.chrom:
                gaps.append(abs(self.second.position - other.second.position))
        elif self.first.chrom == other.second.chrom:
            gaps.append(abs(self.first.position - other.second.position))
            if self.second.chrom == other.first.chrom:
                gaps.append(abs(self.second.position - other.first.position))
        if len(gaps) == 2 and min(gaps) <= consensus_interval:
            return True
        return False

    def is_same_chrom(self) -> bool:
        """Check if fusion have same chromosome on the two breakpoints.

        Return
        ------
        bool
            True if the two partners have the same chromosoms.
        """
        if self._first.is_init() and self._second.is_init():
            if self._first.chrom == self._second.chrom:
                return True
        return False

    def set_region(self, region: region.Region) -> None:
        """Set fusion breakpoint automatically without kwnowledge if the
        first breakpoint is set.

        Return
        ------
        None
        """
        if not self._first.is_init():
            self._first = copy.deepcopy(region)
        else:
            self._second = copy.deepcopy(region)

    def swap_region(self) -> None:
        """Swap the first and second attributes.

        Return
        ------
        None
        """
        tmp = self._first
        self._first = self._second
        self.second = tmp

    # Import Export
    def to_dict(self) -> Dict[Any, Any]:
        """Export object as a dict.

        Return
        ------
        Dict[str, object]
            Object represented by a dictionary.
        """
        return dict(
            first=self._first.to_dict(),
            second=self.second.to_dict(),
            evidence=self.evidence.to_dict(),
            number=self.number,
            is_consensus=self.is_consensus,
            software=self.software,
        )

    @classmethod
    def from_dict(cls, data: Dict[Any, Any]) -> "Fusion":
        """Build object from a dictionary.

        Parameters
        ----------
        data
            A dictionary with all data needed to build a Fusion.

        Return
        ------
        Fusion
            A novel object.
        """
        f = Fusion()
        f.first = region.Region.from_dict(data.get("first", {}))
        f.second = region.Region.from_dict(data.get("second", {}))
        f.evidence = hmn_evidence.Evidence.from_dict(data.get("evidence", {}))
        f.number = data.get("number", 0)
        f.is_consensus = data.get("is_consensus", False)
        f.software = data.get("software", "")
        return f

    # Meta functions
    def __key(self):
        return (self._first, self._second, self._evidence, self._software)

    def __repr__(self) -> str:
        return "%s: %s" % (self.get_name(), str(self.evidence))

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Fusion):
            return self.__key() == other.__key()
        return NotImplemented

    def __lt__(self, other: "Fusion") -> bool:
        return self._evidence < other.evidence

    def __gt__(self, other: "Fusion") -> bool:
        return self._evidence > other.evidence


class FusionEncoder(JSONEncoder):
    """A class to encode a Fusion as a JSON

    Methods
    -------
    default
        Let to encode data
    """

    def default(self, o):
        return o.to_dict()
