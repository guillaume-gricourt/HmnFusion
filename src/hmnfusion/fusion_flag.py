import enum


class FusionFlag(enum.IntEnum):
    """FusionFlag has some class attributes to describe a fusion.

    Class attributes
    ----------
    INTEREST: int=1
        Fusion is interested
    CONSENSUS_PRIMARY: int=2
        Fusion is a consensus at the first level
    CONSENSUS_SECONDARY: int=4
        Fusion is a consensus at the second level
    GENEFUSE: int=8
        Fusion comes from Genefuse
    HMNFUSION: int=16
        Fusion comes from HmnFusion
    LUMPY: int=32
        Fusion comes from Lumpy
    Methods
    -------
    No constructor
    """

    INTEREST = 1
    CONSENSUS_PRIMARY = 2
    CONSENSUS_SECONDARY = 4
    GENEFUSE = 8
    HMNFUSION = 16
    LUMPY = 32

    @classmethod
    def is_interest(cls, flag: int) -> bool:
        return bool(flag & cls.INTEREST)

    @classmethod
    def is_consensus_primary(cls, flag: int) -> bool:
        return bool(flag & cls.CONSENSUS_PRIMARY)

    @classmethod
    def is_consensus_secondary(cls, flag: int) -> bool:
        return bool(flag & cls.CONSENSUS_SECONDARY)

    @classmethod
    def is_consensus(cls, flag: int) -> bool:
        return cls.is_consensus_primary(flag=flag) or cls.is_consensus_secondary(
            flag=flag
        )

    @classmethod
    def is_genefuse(cls, flag: int) -> bool:
        return bool(flag & cls.GENEFUSE)

    @classmethod
    def is_hmnfusion(cls, flag: int) -> bool:
        return bool(flag & cls.HMNFUSION)

    @classmethod
    def is_lumpy(cls, flag: int) -> bool:
        return bool(flag & cls.LUMPY)
