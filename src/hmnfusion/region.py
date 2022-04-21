import operator
import os
import re
import tempfile
from collections import Counter
from typing import Dict, Tuple

import pysam
from hmnfusion import utils


class Region(object):
    """Region let to store genomic location of a one partner of a fusion.

    Attributes
    ----------
    chrom : str
        The chromosome name.
    position: int
        Genomic coordinate.
    start: int
        Start of the interval
    end: int
        End of the interval
    interval: int
        Number of bases
    sequence_reference: str
        Sequence corresponding to the reference
    sequence_sample: str
        Sequence corresponding to the sample
    orientation: str, either "undefined", "left" or "right"
        The orientation of the fusion.

    Methods
    -------
    __init__(chrom="", position=0, start=0, end=0, interval=0, sequence_reference="",
            sequence_sample="", orientation="undefined")
        Construct an object. All arguments are optional.
    """

    def __init__(
        self,
        chrom: str = "",
        position: int = 0,
        start: int = 0,
        end: int = 0,
        interval: int = 0,
        sequence_reference: str = "",
        sequence_sample: str = "",
        orientation: str = "undefined",
    ):
        self._chrom = chrom
        self._position = position
        self._start = start
        self._end = end
        self._interval = interval
        self._sequence_reference = sequence_reference
        self._sequence_sample = sequence_sample
        self._orientation = orientation

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
    def start(self) -> int:
        return self._start

    @start.setter
    def start(self, start: int) -> None:
        self._start = int(start)

    @property
    def end(self) -> int:
        return self._end

    @end.setter
    def end(self, end: int) -> None:
        self._end = int(end)

    @property
    def interval(self) -> int:
        return self._interval

    @interval.setter
    def interval(self, interval: int) -> None:
        self._interval = int(interval)

    @property
    def sequence_reference(self) -> str:
        return self._sequence_reference

    @sequence_reference.setter
    def sequence_reference(self, sequence_reference: str) -> None:
        self._sequence_reference = sequence_reference

    @property
    def sequence_sample(self) -> str:
        return self._sequence_sample

    @sequence_sample.setter
    def sequence_sample(self, sequence_sample: str) -> None:
        self._sequence_sample = sequence_sample

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
    def get_start(self) -> int:
        """Get start attribute. If it has not been initialized, it will be initialize with position and interval attribute.

        Return
        ------
        int
            start attribute
        """
        if self.start > 0:
            return self.start
        if self.position > 0 and self.interval > 0:
            start = self.position - int(self.interval / 2)
            if start < 0:
                start = 0
            self.start = start
        return self.start

    def get_end(self) -> int:
        """Get end attribute. If it has not been initialized, it will be initialize with position and interval attribute.

        Return
        ------
        int
            end attribute
        """
        if self.end > 0:
            return self.end
        if self.position > 0 and self.interval > 0:
            end = self.position + int(self.interval / 2)
            if end < 0:
                end = 0
            self.end = end
        return self.end

    def get_length(self) -> int:
        """Return the length of the region

        Return
        ------
        int
            The length of the region
        """
        if self.interval > 0:
            return self.interval
        if self.end > 0 and self.start >= 0:
            length = self.end - self.start
            if length < 0:
                return 0
            return length
        return 0

    def get_middle(self) -> int:
        """Return half of the interval

        Raises
        ------
        ValueError
            if length of the interval is not an even number

        Return
        ------
        int
            The half of the interval
        """
        if self.get_length() % 2 != 0:
            raise ValueError("Length must be an even number")
        return int(self.get_length() / 2)

    def format(self, fusion: bool = True) -> str:
        """Return a string representation of the object

        Parameters
        ----------
        fusion: bool (default: True)
            if True, format as "<chrom>:<position>"
            else format as "<chrom>:<start>-<end>"

        Return
        ------
        str
            A string representation of the object.
        """
        if fusion:
            return "%s:%s" % (self.chrom, self.position)
        fmt = "%s:%s" % (self.chrom, self.get_start())
        end = self.get_end()
        if end > 0:
            fmt += "-%s" % (end,)
        return fmt

    def is_init(self) -> bool:
        """Check if a Region is initialized.
        Only based on chromosome information

        Return
        ------
        bool
            True if the Region information is set.
        """
        if self.chrom == "":
            return False
        return True

    def set_sequence_reference(self, path: str) -> None:
        """Set sequence_reference from a fasta reference.

        Parameters
        ----------
        path: str
            Path of the fasta file.

        Return
        ------
        bool
            True if the Region information is set.
        """
        rec = pysam.faidx(path, self.format(fusion=False))
        rec = rec.splitlines()
        sequence = "".join(rec[1:])
        sequence = sequence.upper()
        self.sequence_reference = sequence

    def set_sequence_sample(self, path: str) -> None:
        """Set sequence_sample attribute from a bam file. Create a consensus sequence in the interval of the region.

        Parameters
        ----------
        path: str
            Path of the bam file.

        Return
        ------
        None
        """
        tmpfile, stats = Region.filter_bam(path, self)
        self.sequence_sample = Region.create_consensus(tmpfile, self)
        os.remove(tmpfile)

    # Import Export
    def to_dict(self) -> Dict:
        """Export object as a dict.

        Return
        ------
        Dict[str, Union[str, int]]
            Object represented by a dictionary.
        """
        return dict(
            chrom=self.chrom,
            position=self.position,
            start=self.start,
            end=self.end,
            interval=self.interval,
            sequence_reference=self.sequence_reference,
            sequence_sample=self.sequence_sample,
            orientation=self.orientation,
        )

    @classmethod
    def from_dict(cls, data: Dict) -> "Region":
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
        region.chrom = data.get("chrom", "")
        region.position = data.get("position", 0)
        region.start = data.get("start", 0)
        region.end = data.get("end", 0)
        region.interval = data.get("interval", 0)
        region.sequence_reference = data.get("sequence_reference", "")
        region.sequence_sample = data.get("sequence_sample", "")
        region.orientation = data.get("orientation", "")
        return region

    @classmethod
    def from_str(cls, data: str) -> "Region":
        """Build object from a string.

        Parameters
        ----------
        data
            A string formatting as:
            - "<contig>:<position>"
            - "<contig>:<start>-<end>".

        Return
        ------
        Region
            A novel object.
        """
        m = re.search(r"([\d\w]+):(\d+)(-?\d*)", data)
        region = Region()
        if m is not None:
            region.chrom = m.group(1)
            if m.group(3):
                region.start = int(m.group(2))
                region.end = int(m.group(3).replace("-", ""))
            else:
                region.position = int(m.group(2))
        return region

    @classmethod
    def check_region(cls, region: str) -> bool:
        """Verify if a string corresponds to a region

        Parameters
        ----------
        region: str
            A string to parse.

        Return
        ------
        bool:
            True if it corresponds to a region, False otherwise.
        """
        if re.search(r"[\d\w]+:\d+", region, flags=re.I):
            return True
        return False

    @staticmethod
    def filter_bam(
        path: str, region: "Region", soft_clip: int = 4
    ) -> Tuple[str, Dict[str, int]]:
        """Set sequence_sample attribute from a bam file. Create a consensus sequence in the interval of the region.

        Parameters
        ----------
        path: str
            Path of the bam file.
        region: hmnfusion.region.Region
            A region into filter the bam
        soft_clip: int (default: 6)
            Minimum number of base soft clipped to take account

        Return
        ------
        Tuple[str, Dict[Any]]
        """
        alignment = pysam.AlignmentFile(path)
        tmpfile = tempfile.NamedTemporaryFile(delete=False, suffix=".bam")
        stats = dict(total=0, supplementary=0, mate=0, soft_clipped=0)
        with pysam.AlignmentFile(tmpfile.name, "wb", header=alignment.header) as fod:
            for aseg in alignment.fetch(region=region.format(fusion=False)):
                # Filtering.
                if aseg.is_unmapped or aseg.is_duplicate or aseg.is_supplementary:
                    continue
                stats["total"] += 1
                # Count split reads.
                if aseg.has_tag("SA"):
                    fod.write(aseg)
                    stats["supplementary"] += 1
                    continue
                # Count other Chrom.
                if aseg.is_paired:
                    if not aseg.mate_is_unmapped and not aseg.is_unmapped:
                        if aseg.next_reference_id != aseg.reference_id:
                            fod.write(aseg)
                            stats["mate"] += 1
                            continue
                # Count reads clipped.
                for cigar in aseg.cigartuples:
                    if cigar[0] in [4, 5] and cigar[1] >= 6:
                        fod.write(aseg)
                        stats["soft_clipped"] += 1
                        break
        utils.check_bam_index(tmpfile.name)
        return tmpfile.name, stats

    @staticmethod
    def create_consensus(path: str, region: "Region") -> str:
        """Build a sequence consensus into the specified region.

        Parameters
        ----------
        path: str
            Path of the bam file.
        region: hmnfusion.region.Region
            A region into filter the bam.

        Return
        ------
        str
            The consensus sequence.
        """
        alignment = pysam.AlignmentFile(path)
        pos = region.get_start() - 1
        cons = ""
        for pcol in alignment.pileup(
            truncate=True,
            region=region.format(fusion=False),
            max_depth=1000000,
            ignore_overlaps=False,
            ignore_orphans=False,
        ):
            pdepth = pcol.get_num_aligned()  # pileupcolumn.nsegments
            while pos != pcol.reference_pos:
                cons += "N"
                pos += 1
            if pdepth == 0:
                cons += "N"
                pos += 1
                continue
            if pos != pcol.reference_pos:
                raise ValueError(
                    "Shift between reference and pileup : ref %s pileup %s"
                    % (pos, pcol.reference_pos)
                )
            count_base = dict(Counter([x.upper() for x in pcol.get_query_sequences()]))
            if sum(count_base.values()) != pdepth:
                raise ValueError(
                    "Shift between depth: %s and pileup: %s at postion: %s"
                    % (sum(count_base.values()), pdepth, pos)
                )
            cons += max(count_base.items(), key=operator.itemgetter(1))[0]
            pos += 1
        # If not covered at the end of bam
        while pos < region.get_length():
            cons += "N"
            pos += 1
        return cons

    # Meta functions
    def __key(self):
        return (
            self.chrom,
            self.position,
            self.start,
            self.end,
            self.interval,
            self.sequence_reference,
            self.sequence_sample,
            self.orientation,
        )

    def __repr__(self) -> str:
        return "%s %s:%s" % (self.orientation, self.chrom, self.position)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other: object):
        if isinstance(other, Region):
            return self.__key() == other.__key()
        return NotImplemented
