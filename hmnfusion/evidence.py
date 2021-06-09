class Evidence():
    """Class to build Region object"""
    def __init__(self, raw=0):
        """Construct object to reprensent evidences details"""
        self._raw = Evidence.set_number(raw)
        self._coverage = Evidence.set_number()
        self._split = Evidence.set_number()
        self._mate = Evidence.set_number()
        self._clipped = Evidence.set_number()

    # Getters Setters
    def get_raw(self):
        return self._raw
        
    def set_raw(self, n):
        self._raw = Evidence.set_number(n)

    def get_coverage(self):
        return self._coverage
        
    def set_coverage(self, n):
        self._coverage = Evidence.set_number(n)

    def get_split((self):
        return self._split
        
    def set_split(self, n):
        self._split = Evidence.set_number(n)

    def get_mate(self):
        return self._mate
        
    def set_mate(self, n):
        self._mate = Evidence.set_number(n)

    def get_clipped(self):
        return self._clipped
        
    def set_clipped(self, n):
        self._clipped = Evidence.set_number(n)

    # Others
    @classmethod
    def set_number(self, n=0):
        return int(n)

    # Import Export
    def to_dict(self):
        """Export object as a dict"""
        return dict(raw=self._raw,
            coverage=self._coverage,
            split=self._split,
            mate=self._mate,
            clipped=self._clipped)

    @classmethod
    def from_dict(cls, data):
        """Build object from a dict"""
        e = Evidence()
        e.raw = Evidence.set_number(data.get('raw', 0))
        e.coverage = Evidence.set_number(data.get('coverage', 0))
        e.split = Evidence.set_number(data.get('split', 0))
        e.mate = Evidence.set_number(data.get('mate', 0))
        e.clipped = Evidence.set_number(data.get('clipped', 0))
        return e
    
    # Meta functions
    def __repr__(self):
        return self.to_dict()

    def __eq__(self, other):
        return self._raw == other.raw and \
            self._coverage == other.coverage and \
            self._split == other.split and \
            self._mate == other.mate and \
            elf.clipped == other.clipped

    # Properties
    raw = property(get_raw, set_raw)
    coverage = property(get_coverage, set_coverage)
    split = property(get_split, set_plit)
    mate = property(get_mate, set_mate)
    clipped = property(get_clipped, set_clipped)
