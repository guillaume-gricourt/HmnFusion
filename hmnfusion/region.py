class Region():
    """Class to build Region object"""
    def __init__(self, chrom='', position=0, orientation='undefined'):
        """Construct object with a chromosome, a position and an orientation"""
        self.set_chrom(chrom)
        self.set_position(position)
        self.set_orientation(orientation)

    # Getters Setters
    def get_chrom(self):
        return self._chrom
        
    def set_chrom(self, chrom):
        self._chrom = chrom
        
    def get_position(self):
        return self._position
        
    def set_position(self, position):
        self._position = int(position)

    def get_orientation(self):
        return self._orientation
    
    def set_orientation(self, orientation):
        if orientation in ['left', 'right']:
            self._orientation = orientation
        else:
            self._orientation = 'undefined'

    # Others
    def is_init(self):
        """Chect if a Region is initialize. Only based on chromosome information"""
        if self._chrom == '':
            return False
        return True

    # Import Export
    def to_dict(self):
        """Export object as a dict"""
        return dict(orientation=self._orientation,
            chrom=self._chrom,
            position=self._position)

    @classmethod
    def from_dict(cls, data):
        """Build object from a dict"""
        region = Region()
        region.orientation = data.get('orientation', '')
        region.chrom = data.get('chrom', '')
        region.position = data.get('position', 0)
        return region
    
    # Meta functions
    def __repr__(self):
        return '%s %s:%s'%(self._orientation, self._chrom, self._position)

    def __eq__(self, other):
        return self._chrom == other.chrom and \
            self._position == other.position and \
            self._orientation == other.orientation

    # Properties
    chrom = property(get_chrom, set_chrom)
    position = property(get_position, set_position)
    orientation = property(get_orientation, set_orientation)
