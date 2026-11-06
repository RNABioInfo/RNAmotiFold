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


class AvailableBinaries(StrEnum):
    RNAMOTIFOLD = "RNAmotiFold"
    RNAMOSH = "RNAmoSh"
    RNAMOTICES = "RNAmotiCes"
    RNAMOTIALIGN = "RNAmotiAlign"
    RNAMOTIFOLD_MOTIFLP = "RNAmotiFoldMotmicro"
    RNAMOSH_MOTIFLP = "RNAmoShMotmicro"
    RNAMOTICES_MOTIFLP = "RNAmotiCesMotmicro"
    RNAMOTIFOLD_PFC = "RNAmotiFold_pfc"
    RNAMOSH_PFC = "RNAmoSh_pfc"
    RNAMOTICES_PFC = "RNAmotiCes_pfc"
    RNAMOTIFOLD_MOTIFLP_PFC = "RNAmotiFold_motmacro_pfc"
    RNAMOSH_MOTIFLP_PFC = "RNAmoSh_motmacro_pfc"
    RNAMOTICES_MOTIFLP_PFC = "RNAmotiCes_motmacro_pfc"
    RNAMOTIFOLD_SUBOPT = "RNAmotiFold_subopt"
    RNAMOSH_SUBOPT = "RNAmoSh_subopt"
    RNAMOTICES_SUBOPT = "RNAmotiCes_subopt"
    RNAMOTIFOLD_MOTIFLP_SUBOPT = "RNAmotiFold_motmacro_subopt"
    RNAMOSH_MOTIFLP_SUBOPT = "RNAmoSh_motmacro_subopt"
    RNAMOTICES_MOTIFLP_SUBOPT = "RNAmotiCes_motmacro_subopt"


def GetProcessType(Algorithm: AlgorithmType) -> ProcessType:
    match Algorithm:
        case AlgorithmType.RNAMOTIALIGN:
            return ProcessType.ALIGNMENT
        case (
            AlgorithmType.RNAMOTIFOLD
            | AlgorithmType.RNAMOTICES
            | AlgorithmType.RNAMOSH
        ):
            return ProcessType.SINGLE
