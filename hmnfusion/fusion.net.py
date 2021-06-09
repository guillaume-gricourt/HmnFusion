import networkx nx

from networkx.readwrite import json_graph
from .region import Region
    
class Fusion():
    """Class Fusion. Provide attributes and methods to manipulate Fusion items"""
    name_consensus = 'CONS'

    def __init__(self):
        """Construct Fusion object with a software name (optional)"""
        self._graph = nx.Graph()
        self._first = Region()
        self._second = Region()
        self._evidence = 0
        self._evidence_details = {}
        self._depth = 0
        self._ident = ''
        self._isConsensus = False
        self._buildFrom = set()
        self._software = set()
        if software:
            self.set_software(software)
    
    # Getters Setters
    def add_node(self):
        return self._first
        
    # Others
    def is_near(self, other, consensus_interval):
        """Check if fusion are near from an other with a distance parameter"""
        gaps = []
        if self._first.chrom == other.first.chrom:
            gaps.append(abs(self._first.position-other.first.position))
            if self._second.chrom == other.second.chrom:
                gaps.append(abs(self._second.position-other.second.position))
        elif self._first.chrom == other.second.chrom:
            gaps.append(abs(self._first.position-other.second.position))
            if self._second.chrom == other.first.chrom:
                gaps.append(abs(self._second.position-other.first.position))
        if len(gaps) == 2 and min(gaps) <= consensus_interval:
            return True
        return False

    def remove_fom_name_cons(self):
        """Remove consensus name in buildFrom attribute"""    
        self._buildFrom = set([ x for x in self._buildFrom if not x.startswith(Fusion.name_consensus)])

    def is_same_chrom(self):
        """Check if fusion have same chromosome on the two breakpoints"""
        if self._first.is_init() and self._second.is_init():
            if self._first.chrom == self._second.chrom:
                return True
        return False

    def set_region(self, region):
        """Set fusion breakpoint automatically without kwnowledge if the first breakpoint is set"""
        if not self._first.is_init():
            self._first = region
        else:
            self._second = region

    def swap_region(self):
        """Swap the breakpoints"""    
        tmp = self._first
        self._first = self._second
        self.second = tmp

    # Import Export
    def to_dict(self):
        """Export Fusion as dict"""
        return dict(software=list(self._software), 
                first=self._first.to_dict(), 
                second=self._second.to_dict(), 
                evidence=self._evidence,
                evidence_details=self._evidence_details,
                ident=self._ident, 
                buildFrom=list(self._buildFrom),
                isConsensus=self._isConsensus,
                depth=self._depth)

        self._first = Region()
        self._second = Region()
        self._evidence = 0
        self._evidence_details = {}
        self._depth = 0
        self._ident = ''
        self._isConsensus = False
        self._buildFrom = set()
        self._software = set()
        if software:
            self.set_software(software)

    @classmethod
    def from_dict(cls, data):
        """Construct a Fusion object from a dict"""    
        fusion = Fusion()
        fusion.first = Region.from_dict(data.get('first', {}))
        fusion.second = Region.from_dict(data.get('second', {}))
        fusion.evidence = data.get('evidence', 0)
        fusion.evidence_details = data.get('evidence_details', {})
        fusion.depth = data.get('depth', 0)
        fusion.ident = data.get('ident', '')
        fusion.isConsensus = data.get('isConsensus', False)
        fusion.buildFrom = set(data.get('buildFrom', []))
        fusion.software = set(data.get('software', []))
        return fusion

    # Meta functions
    def __repr__(self):
        return 'From %s %s %s (%s) %s %s Evidence %s Depth %s'%(','.join(self._software), self._ident, ','.join(self._buildFrom), self._evidence, self._first, self._second, self._evidence, self._depth)
 
    def __eq__(self, other):
        return self._first == other.first and \
            self._second == other.second and \
            self._evidence == other.evidence and \
            self._evidence_details == other.evidence_details and \
            self._depth == other.depth and \
            self._isConsensus == other.isConsensus and \
            self._ident == other.ident and \
            len(self._buildFrom.union(other.buildFrom)) == len(self._buildFrom) and \
            len(self._buildFrom.union(other.buildFrom)) == len(other.buildFrom) and \
            len(self._software.union(other.software)) == len(self._software) and \
            len(self._software.union(other.software)) == len(other.software)

    def __lt__(self, other):
        return self._evidence < other.evidence

    def __gt__(self, other):
        return self._evidence > other.evidence

    # Properties
    first = property(get_first, set_first)
    second = property(get_second, set_second)
    evidence = property(get_evidence, set_evidence)
    evidence_details = property(get_evidence_details, set_evidence_details)
    depth = property(get_depth, set_depth)
    ident = property(get_ident, set_ident)
    buildFrom = property(get_build_from, set_build_from)
    isConsensus = property(get_is_consensus, set_is_consensus)
    software = property(get_software, set_software)
