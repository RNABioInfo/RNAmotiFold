import logging
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from RNAmotiFold.input.parameters import ScriptParameters
from RNAmotiFold.Misc.EnumTypes import (
    AlgorithmType,
    AllowLonelyBP,
    MotifOrientation,
    MotifSource,
    ShapeLevel,
)
from RNAmotiFold.Misc.ValueObjects import (
    SubOptEnergy,
    fraction,
    limited_temperature,
    pos_int,
)

from RNAmotiFold import RNAMOTIFOLD_ROOT_DIR

logger = logging.getLogger(__name__)

# A Python class for making Bellman's GAP more convenient to use
# Just create a class instances, feed it with the call arguments you need
# and it'll create a call from the arguments classified as runtime arguments.
# These are: motif_source, motif_orientation, kvalue, hishape_mode, shape_level, energy


def _set_algorithm_path(algorithm_binary: str) -> Path:
    path = RNAMOTIFOLD_ROOT_DIR / "bin" / algorithm_binary
    logger.debug(f"Set algorithm path as {path}")
    return path


def _set_algorithm_binary(
    algorithm: AlgorithmType,
    subopt: bool,
    pfc: bool,
    basepairs: AllowLonelyBP,
) -> str:
    if algorithm == AlgorithmType.RNAMOTIALIGN:
        if subopt or pfc:
            raise NotImplementedError(
                "RNAmotiAlign is not available with subopt or partition function implementations"
            )
        return "RNAmotiAlign"
    match (basepairs, subopt, pfc):
        case (2, False, False):
            alg = algorithm.value + "Motmicro"
        case (2, True, False):
            alg = algorithm.value + "_motmacro_subopt"
        case (2, False, True):
            alg = algorithm.value + "_motmacro_pfc"
        case (1, False, False) | (0, False, False):
            alg = algorithm.value
        case (1, True, False) | (0, True, False):
            alg = algorithm.value + "_subopt"
        case (1, False, True) | (0, False, True):
            alg = algorithm.value + "_pfc"
        case _:
            raise ValueError(
                "The algorithm you specified does not exist, please revisit your arguments."
            )
    logger.debug(f"Set algorithm binary as {alg}")
    return alg


@dataclass
class MainCall:
    algorithm_path: Path
    runtime_dict: dict[str, str]

    def to_list(self) -> list[str]:
        arguments = [str(self.algorithm_path)]
        for k, v in self.runtime_dict.items():
            if v is not None:
                arguments.extend([k, str(v)])
        return arguments


# I should really make this function just take the parameters and then have a wrapper
# For Calling it from ScriptParams or Prediction Objects


def make_main_call(
    motif_source: MotifSource,
    motif_orientation: MotifOrientation,
    temperature: limited_temperature,
    version: str,
    pfc: bool,
    subopt: bool,
    subopt_opts: SubOptEnergy,
    low_prob_filter: fraction,
    kvalue: pos_int,
    algorithm: AlgorithmType,
    shape_level: ShapeLevel,
    motif_weight: Decimal,
    motif_fraction: fraction,
    basepairs: AllowLonelyBP,
) -> MainCall:
    """Automatic call setter. The function checks the set algorithm and builds a call string based on it."""

    runtime_dictionary: dict[str, str] = {
        "-Q": str(motif_source.value),
        "-b": str(motif_orientation.value),
        "-t": str(temperature),
        "-v": str(version),
    }
    # Only one of these will be activate at a time
    if subopt:
        runtime_dictionary[subopt_opts.commandline_parameter] = str(
            subopt_opts.subopt_energy
        )
    elif pfc:
        runtime_dictionary["-F"] = str(low_prob_filter)
    else:
        runtime_dictionary["-k"] = str(kvalue)

    if algorithm == AlgorithmType.RNAMOSH:
        runtime_dictionary["-q"] = str(shape_level)
    if algorithm == AlgorithmType.RNAMOTIALIGN:
        runtime_dictionary["-W"] = str(motif_weight)
        runtime_dictionary["-D"] = str(motif_fraction)
    if basepairs in [AllowLonelyBP.NO_LP, AllowLonelyBP.ALL_LP]:
        runtime_dictionary["-u"] = str(basepairs.value)
    algorithm_binary = _set_algorithm_binary(
        algorithm, subopt, pfc, basepairs
    )
    mc = MainCall(_set_algorithm_path(algorithm_binary), runtime_dictionary)
    logger.debug(f"Set call as {' '.join(mc.to_list())}")
    return mc


def main_call_from_script_params(rt_args: ScriptParameters) -> MainCall:
    """Automatic call setter from ScriptParameters. The function checks the set algorithm and builds a call string based on it."""
    return make_main_call(
        rt_args.motif_source,
        rt_args.motif_orientation,
        rt_args.temperature,
        rt_args.version,
        rt_args.pfc,
        rt_args.subopt,
        rt_args.subopt_opts,
        rt_args.low_prob_filter,
        rt_args.kvalue,
        rt_args.algorithm,
        rt_args.shape_level,
        rt_args.motif_weight,
        rt_args.motif_fraction,
        rt_args.basepairs,
    )
