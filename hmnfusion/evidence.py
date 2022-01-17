class Evidence(object):
    """Class to build Region object"""
    def __init__(self, raw=0):
        """Construct object to reprensent evidences details"""
        self._raw = Evidence.set_number(raw)
        self._split = Evidence.set_number()
        self._mate = Evidence.set_number()
        self._clipped = Evidence.set_number()
        self._depth = Evidence.set_number()

    # Getters Setters
    @property
    def raw(self):
        return self._raw

    @raw.setter
    def raw(self, n):
        self._raw = Evidence.set_number(n)

    @property
    def split(self):
        return self._split

    @split.setter
    def split(self, n):
        self._split = Evidence.set_number(n)

    @property
    def mate(self):
        return self._mate

    @mate.setter
    def mate(self, n):
        self._mate = Evidence.set_number(n)

    @property
    def clipped(self):
        return self._clipped

    @clipped.setter
    def clipped(self, n):
        self._clipped = Evidence.set_number(n)

    @property
    def depth(self):
        return self._depth

    @depth.setter
    def depth(self, n):
        self._depth = Evidence.set_number(n)

    # Others
    @classmethod
    def set_number(self, n=0):
        return abs(int(float(n)))

    def get_sum(self):
        return self._split + self._mate + self._clipped

    def get_vaf(self):
        vaf = 0
        if self._depth > 0:
            vaf = (self.get_sum() / self._depth) * 100
        return '{:.2f}'.format(vaf).replace('.', ',')

    def get_max_count(self):
        return max([self._raw, self.get_sum()])

    # Import Export
    def to_dict(self):
        """Export object as a dict"""
        return dict(
            raw=self._raw,
            split=self._split,
            mate=self._mate,
            clipped=self._clipped,
            depth=self._depth
        )

    @classmethod
    def from_dict(cls, data):
        """Build object from a dict"""
        e = Evidence()
        e.raw = Evidence.set_number(data.get('raw', 0))
        e.split = Evidence.set_number(data.get('split', 0))
        e.mate = Evidence.set_number(data.get('mate', 0))
        e.clipped = Evidence.set_number(data.get('clipped', 0))
        e.depth = Evidence.set_number(data.get('depth', 0))
        return e

    # Meta functions
    def __key(self):
        return (
            self._raw,
            self._split,
            self._mate,
            self._clipped,
            self._depth
        )

    def __repr__(self):
        return 'Ev ' + ' '.join(str(x) for x in self.__key())

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Evidence):
            return self.__key() == other.__key()
        return NotImplemented

    def __lt__(self, other):
        return self.get_max_count() < other.get_max_count()

    def __gt__(self, other):
        return self.get_max_count() > other.get_max_count()
