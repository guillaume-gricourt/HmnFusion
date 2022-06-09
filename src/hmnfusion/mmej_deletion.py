import hashlib
from typing import List

import pandas as pd
import pysam
from hmnfusion import utils
from natsort import natsorted


class Conclude(utils.EnumNoValue):
    """Conclude is an Enum class to represent conclusion about the mmej sequence found"""

    UNINITIALIZED = ""
    AMBIGUOUS = "alignment ambiguous"
    UNCLEAR = "no clear signature"
    VALID = "mmej signature"


class Value(object):
    """Value store a record corresponding to one Point of a MmejDeletion

    Attributes
    ----------
    id: str
        The id of the object
    contig: str
        Contig id
    start: int
        Start of the deletion
    deletion: str
        Sequence corresponding to the deletion found in the sample
    sequence: str
        Sequence corresponding to the mmej

    Methods
    -------
    __init__(id="", position="", start=0, deletion="", sequence="")
        Construct an object. All arguments are optional.
    """

    def __init__(
        self,
        id: str = "",
        contig: str = "",
        start: int = 0,
        deletion: str = "",
        sequence: str = "",
    ) -> None:
        self.id = id
        self.contig = contig
        self.start = int(start)
        self.deletion = deletion
        self.sequence = sequence

    def get_conclusion(self) -> Conclude:
        """Build conclusion based on the sequence and deletion attributes

        Return
        ------
        Conclude
            The appropriate conclusion
        """
        len_deletion = len(self.deletion)
        len_mh = len(self.sequence)
        # if len_deletion <= 1 or len_mh <= 1:
        # pass
        if len_deletion == len_mh:
            return Conclude.AMBIGUOUS
        elif len_deletion > 1 and len_mh < 5:
            return Conclude.UNCLEAR
        elif len_deletion > 1 and len_mh >= 5:
            return Conclude.VALID
        return Conclude.UNINITIALIZED

    def set_sequence(self, path: str) -> None:
        """Set the sequence attribute

        Parameters
        ----------
        path: str
            Path of the fasta reference to extract the information

        Return
        ------
        None
        """
        len_deletion = len(self.deletion)
        faidx = pysam.faidx(path, self.to_region(), split_lines=True)
        seq = "".join(faidx[1:]).upper()

        left = seq[:len_deletion]
        right = seq[len_deletion:]

        assert self.deletion == left

        mh_len = 0
        for i in range(len_deletion + 1):
            motif = left[:i]
            if right.startswith(motif):
                mh_len = i
                self.sequence = motif
            if i > mh_len:
                break

    @classmethod
    def from_record(cls, record: pysam.VariantRecord) -> "Value":
        """Build a Value object from a pysam.VariantRecord

        Parameters
        ----------
        record: pysam.VariantRecord
            An object to transform

        Return
        ------
        Value
            The novel object
        """
        record_name = "%s %s %s %s" % (
            record.contig,
            record.pos,
            record.ref,
            record.alts[0],
        )
        record_id = hashlib.md5(record_name.encode("utf8")).hexdigest()
        return Value(
            id=record_id,
            contig=record.contig,
            start=record.pos + 1,
            deletion=record.ref[1:].upper(),
        )

    def to_dataframe(self) -> pd.DataFrame:
        """Convert the value object to a dataframe

        Return
        ------
        pd.DataFrame
            A dataframe with five columns:
                - contig
                - start
                - deletion
                - sequence
                - conclusion
            and one index:
                - id of the object
        """
        return pd.DataFrame(
            {
                "contig": self.contig,
                "start": self.start,
                "deletion": self.deletion,
                "sequence": self.sequence,
                "conclusion": self.get_conclusion().value,
            },
            index=[self.id],
        )

    def to_region(self) -> str:
        """Convert a Value object into a region string

        Return
        ------
        str
            The formatting region
        """
        return "%s:%s-%s" % (
            self.contig,
            self.start,
            self.start + (2 * len(self.deletion)),
        )

    # Meta functions
    def __key(self):
        return (
            self.id,
            self.contig,
            self.start,
            self.deletion,
            self.sequence,
        )

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other: object):
        if isinstance(other, Value):
            return self.__key() == other.__key()
        return NotImplemented


class MmejDeletion(object):
    """MmejDeletion is a class to deal for deletion mmej

    Attributes
    ----------
    name: str
        The name of the sample
    values: List[Value]
        List of records

    Methods
    -------
    __init__(name, values)
        Construct an object. Name is required.
    """

    # Be careful, if you use an empty list as default args
    # it would be shared between objects
    def __init__(
        self,
        name: str,
        values: List[Value],
    ) -> None:
        self.name = name
        self.values = values

    @property
    def empty(self) -> bool:
        if len(self.values) > 0:
            return False
        return True

    @classmethod
    def build_empty_dataframe(cls, name: str) -> pd.DataFrame:
        """Construct an empty dataframe

        Parameters
        ----------
        name: str
            Name of one column to add to the header

        Return
        ------
        pd.DataFrame
            A dataframe with five columns:
                - contig
                - start
                - deletion
                - sequence
                - conclusion
            and one index:
                - id of the object
        """
        return pd.DataFrame(
            columns=["contig", "start", "deletion", "sequence", "conclusion", name]
        )

    def get_value_ids(self) -> List[str]:
        """Return all ids of record

        Return
        ------
        List[str]
            All ids in the records
        """
        return [x.id for x in self.values]

    def set_value_sequence(self, path: str) -> None:
        """Apply set_sequence for each value

        Parameters
        ----------
        path: str

        Return
        ------
        None

        See Also
        --------
        Value().set_sequence()
        """
        for value in self.values:
            value.set_sequence(path=path)

    @classmethod
    def from_vcf(cls, path: str) -> List["MmejDeletion"]:
        """Load an object MmejDeletion from a VCF file

        Parameters
        ----------
        path: str
            Path of the VCF file

        Return
        ------
        List[MmejDeletion]
            Return a list of objects, because a VCF file should contained one or more samples
        """
        # Init
        data = {}
        vcf_in = pysam.VariantFile(path)
        # Populate samples
        for sample in vcf_in.header.samples:
            data[sample] = cls(name=sample, values=[])

        # Loop over recors
        for record in vcf_in.fetch():
            # Check if variant is a deletion.
            if len(record.ref) > len(record.alts[0]):
                value = Value.from_record(record)
                for sample, variant_record_sample in record.samples.items():
                    initialized = set([0])
                    for field in variant_record_sample.values():
                        if field not in [(None,), None]:
                            initialized.add(1)
                    if (
                        len(initialized) > 1
                        and value.id not in data[sample].get_value_ids()
                    ):
                        data[sample].values.append(value)
        return list(data.values())

    @classmethod
    def to_dataframe(cls, mmej_deletions: List["MmejDeletion"]) -> pd.DataFrame:
        """Export a List of this object to a dataframe

        Parameters
        ----------
        mmej_deletions: List[MmejDeletion]
            A list of object to represent

        Return
        ------
        pd.DataFrame
            A dataframe
        """
        df = pd.DataFrame()
        # Samples
        mmej_deletions = natsorted(mmej_deletions, key=lambda x: x.name)
        # Populate index
        for mmej_del in mmej_deletions:
            if mmej_del.empty:
                df = pd.concat([df, cls.build_empty_dataframe(name=mmej_del.name)])
                continue
            for value in mmej_del.values:
                df = pd.concat([df, value.to_dataframe()])
        if df.empty:
            return df
        df = df[~df.index.duplicated(keep="first")]
        # Sort values.
        idx, *_ = zip(
            *natsorted(
                zip(df.index, df.contig, df.start, df.deletion),
                key=lambda x: (x[1], x[2], x[3]),
            )
        )
        df = df.loc[list(idx)]
        # Populate by samples
        for mmej_del in mmej_deletions:
            for value in mmej_del.values:
                df.at[value.id, mmej_del.name] = "o"
        df = df.fillna(pd.NA)
        return df

    @classmethod
    def to_excel(
        cls,
        path: str,
        mmej_deletions: List["MmejDeletion"],
        sheet_name: str = "mmej_deletion",
    ) -> None:
        """Write MmejDeletion objects to an excel file

        Parameters
        ----------
        path: str
            Path of an output file (xlsx, excel format)
        mmej_deletions: List[MmejDeletion]
            A list of the MmejDeletion to write

        Return
        ------
        None
        """
        df = cls.to_dataframe(mmej_deletions)
        writer = pd.ExcelWriter(path)

        if df.empty:
            for header in df.columns:
                df.at["no deletion found", header] = pd.NA
            df.to_excel(writer, sheet_name=sheet_name)
            writer.save()
            return

        # Write output.
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        ws = writer.sheets[sheet_name]  # pull worksheet object
        # Adjust width.
        utils.adjust_dim_worksheet(ws)
        writer.save()

    # Meta functions
    def __key(self):
        return (
            self.name,
            self.values,
        )

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other: object):
        if isinstance(other, MmejDeletion):
            return self.__key() == other.__key()
        return NotImplemented
