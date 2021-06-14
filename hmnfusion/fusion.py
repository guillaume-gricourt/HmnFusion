import copy

from json import JSONEncoder

from . import (evidence, region, _version)
    
class Fusion():
    """Class Fusion. Provide attributes and methods to manipulate Fusion items"""
    def __init__(self, software=_version.__app_name__):
        """Construct Fusion object with a software name (optional)"""
        self._first = region.Region()
        self._second = region.Region()
        self._evidence = evidence.Evidence()
        self._number = 0
        self._is_consensus = False
        self._software = software
    
    # Getters Setters
    @property
    def first(self):
        return self._first
    @first.setter    
    def first(self, first):
        self._first = copy.deepcopy(first)
    @property
    def second(self):
        return self._second
    @second.setter    
    def second(self, second):
        self._second = copy.deepcopy(second)
    @property
    def evidence(self):
        return self._evidence
    @evidence.setter
    def evidence(self, e):
        self._evidence = copy.deepcopy(e)
    @property
    def number(self):
        return self._number
    @number.setter
    def number(self, n):
        self._number = n
    @property
    def is_consensus(self):
        return self._is_consensus
    @is_consensus.setter
    def is_consensus(self, b):
        self._is_consensus = b
    @property
    def software(self):
        return self._software
    @software.setter
    def software(self, software):
        self._software = software

    # Others
    def update(self, other):
        self.first = other.first
        self.second = other.second
        self.evidence = other.evidence
        self.software = other.software

    def get_name(self):
        n = self.software[:3].upper()
        if self.is_consensus:
            n = _version.__app_name__[:3].upper()
        return '%s_%s'%(n ,self._number)

    def is_near(self, other, consensus_interval):
        """Check if fusion are near from an other with a distance parameter"""
        gaps = []
        if self.first.chrom == other.first.chrom:
            gaps.append(abs(self.first.position-other.first.position))
            if self.second.chrom == other.second.chrom:
                gaps.append(abs(self.second.position-other.second.position))
        elif self.first.chrom == other.second.chrom:
            gaps.append(abs(self.first.position-other.second.position))
            if self.second.chrom == other.first.chrom:
                gaps.append(abs(self.second.position-other.first.position))
        if len(gaps) == 2 and min(gaps) <= consensus_interval:
            return True
        return False

    def is_same_chrom(self):
        """Check if fusion have same chromosome on the two breakpoints"""
        if self._first.is_init() and self._second.is_init():
            if self._first.chrom == self._second.chrom:
                return True
        return False

    def set_region(self, region):
        """Set fusion breakpoint automatically without kwnowledge if the first breakpoint is set"""
        if not self._first.is_init():
            self._first = copy.deepcopy(region)
        else:
            self._second = copy.deepcopy(region)

    def swap_region(self):
        """Swap the breakpoints"""    
        tmp = self._first
        self._first = self._second
        self.second = tmp

    # Import Export
    def to_dict(self):
        """Export Fusion as dict"""
        return dict(first=self._first.to_dict(), 
                second=self._second.to_dict(), 
                evidence=self._evidence.to_dict(),
                number=self._number,
                is_consensus=self._is_consensus,
                software=self._software)

    @classmethod
    def from_dict(cls, data):
        """Construct a Fusion object from a dict"""    
        f = Fusion()
        f.first = region.Region.from_dict(data.get('first', {}))
        f.second = region.Region.from_dict(data.get('second', {}))
        f.evidence = evidence.Evidence.from_dict(data.get('evidence', {}))
        f.number = data.get('number', 0)
        f.is_consensus = data.get('is_consensus', False)
        f.software = data.get('software', '')
        return f

    # Meta functions
    def __key(self):
        return (self._first, self._second, self._evidence, self._software)

    def __repr__(self):
        #return 'From %s %s %s %s %s %s %s %s'%(self._software, str(self._first), str(self._second), str(self._evidence), self._number, self._level, self.get_build_from(), self._is_interest)
        return '%s: %s'%(self.get_name(), str(self.evidence))

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Fusion):
            return self.__key() == other.__key()
        return NotImplemented

    def __lt__(self, other):
        return self._evidence < other.evidence

    def __gt__(self, other):
        return self._evidence > other.evidence

class FusionEncoder(JSONEncoder):
    def default(self, o):
        return o.to_dict()
