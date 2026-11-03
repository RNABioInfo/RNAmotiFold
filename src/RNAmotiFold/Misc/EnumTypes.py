from enum import IntEnum, StrEnum


class ProcessType(StrEnum):
    ALIGNMENT = "ali"
    SINGLE = "single"


class AlgorithmType(StrEnum):
    RNAMOTIFOLD = "RNAmotiFold"
    RNAMOTICES = "RNAmotiCes"
    RNAMOSH = "RNAmoSh"
    RNAMOTIALIGN = "RNAmotiAlign"


class PredictionType(StrEnum):
    MFE = "mfe"
    PFC = "pfc"
    ALI = "ali"


class MotifSource(IntEnum):
    BGSU = 1
    RFAM = 2
    ALL = 3


class MotifOrientation(IntEnum):
    FORWARD = 1
    REVERSE = 2
    BOTH = 3


class ShapeLevel(IntEnum):
    LEVEL1 = 1
    LEVEL2 = 2
    LEVEL3 = 3
    LEVEL4 = 4
    LEVEL5 = 5


class AllowLonelyBP(IntEnum):
    NO_LP = 0
    ALL_LP = 1
    MOTIF_LP = 2


class LoopType(StrEnum):
    HAIRPIN = "hairpin"
    INTERNAL = "internal"
    BULGE = "bulge"
    ALL = "all"


class MotifOrigin(StrEnum):
    DEFAULT = "default"
    CUSTOM = "custom"


def GetProcessType(Algorithm:AlgorithmType) -> ProcessType:
    match Algorithm:
        case AlgorithmType.RNAMOTIALIGN:
            return ProcessType.ALIGNMENT
        case AlgorithmType.RNAMOTIFOLD | AlgorithmType.RNAMOTICES | AlgorithmType.RNAMOSH:
            return ProcessType.SINGLE   