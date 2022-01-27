class Region():
    """Class to build Region object"""

    def __init__(self, chrom="", position=0, orientation="undefined"):
        """Construct object with a chromosome, a position and an orientation"""
        self.chrom = chrom
        self.position = position
        self.orientation = orientation

    # Getters Setters
    @property
    def chrom(self):
        return self._chrom

    @chrom.setter
    def chrom(self, chrom):
        self._chrom = chrom

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, position):
        self._position = int(position)

    @property
    def orientation(self):
        return self._orientation

    @orientation.setter
    def orientation(self, orientation):
        if orientation in ["left", "right"]:
            self._orientation = orientation
        else:
            self._orientation = "undefined"

    # Others
    def is_init(self):
        """Chect if a Region is initialize.
        Only based on chromosome information
        """
        if self._chrom == "":
            return False
        return True

    # Import Export
    def to_dict(self):
        """Export object as a dict"""
        return dict(
            orientation=self._orientation, chrom=self._chrom, position=self._position
        )

    @classmethod
    def from_dict(cls, data):
        """Build object from a dict"""
        region = Region()
        region.orientation = data.get("orientation", "")
        region.chrom = data.get("chrom", "")
        region.position = data.get("position", 0)
        return region

    # Meta functions
    def __key(self):
        return (self._orientation, self._chrom, self._position)

    def __repr__(self):
        return "%s %s:%s" % (self._orientation, self._chrom, self._position)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, Region):
            return self.__key() == other.__key()
        return NotImplemented
