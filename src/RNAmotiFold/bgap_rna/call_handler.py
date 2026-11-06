import logging
from pathlib import Path
from typing import Self

from numpy import float128
from RNAmotiFold.input.parameters import ScriptParameters
from RNAmotiFold.Misc.EnumTypes import (
    AlgorithmType,
    AllowLonelyBP,
    MotifOrientation,
    MotifSource,
    ShapeLevel,
)
from RNAmotiFold.Misc.Helpers import interval_check

from RNAmotiFold import RNAMOTIFOLD_ROOT_DIR

logger = logging.getLogger(__name__)

# A Python class for making Bellman's GAP more convenient to use
# Just create a class instances, feed it with the call arguments you need
# and it'll create a call from the arguments classified as runtime arguments.
# These are: motif_source, motif_orientation, kvalue, hishape_mode, shape_level, energy


class CallHandler:
    """Main class for creating the call to the RNAmotiFold, RNAmotiCes and RNAmotiAlign binaries. Creates the correct binary name, algorithm path and organizes commandline arguments"""

    def __repr__(self) -> str:
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(len(v)):
            together.append(f"{k[i]}={v[i]!r}")
        return f"{classname}({', '.join(together)})"

    @classmethod
    def from_script_parameters(cls, params: ScriptParameters) -> Self:
        obj = cls(
            alg=params.algorithm,
            motif_source=params.motif_source,
            motif_orientation=params.motif_orientation,
            kvalue=params.kvalue,
            shape_level=params.shape_level,
            energy=params.energy,
            pfc=params.pfc,
            low_prob_filter=params.low_prob_filter,
            temperature=params.temperature,
            energy_percent=params.energy_percent,
            allowLonelyBasepairs=params.basepairs,
            subopt=params.subopt,
            weight=params.motif_weight,
            fraction=params.motif_fraction,
            version=params.version,
        )
        logger.debug(f"Created bgap obj {obj!s} from {vars(params)}")
        return obj

    def __init__(
        self,
        alg: AlgorithmType = AlgorithmType.RNAMOTIFOLD,
        motif_source: MotifSource = MotifSource.BGSU,
        motif_orientation: MotifOrientation = MotifOrientation.FORWARD,
        kvalue: int = 5,
        shape_level: ShapeLevel = ShapeLevel.LEVEL3,
        energy: float | None = None,
        temperature: float = 37.0,
        energy_percent: float = 5.0,
        allowLonelyBasepairs: AllowLonelyBP = AllowLonelyBP.NO_LP,
        subopt: bool = False,
        pfc: bool = False,
        low_prob_filter: float = 0.000001,
        weight: float = 1.0,
        fraction: float = 0.7,
        version: str = "4_10",
    ) -> None:
        self.subopt = subopt
        self.pfc = pfc
        self.low_probability_filter = low_prob_filter
        self.motif_source = motif_source
        self.motif_orientation = motif_orientation
        self.kvalue = kvalue
        self.shape_level = shape_level
        self.absolute_energy = energy
        self.temperature = temperature
        self.energy_percent = energy_percent
        self.allowLonelyBasepairs = allowLonelyBasepairs
        self.algorithm = alg  # Set algorithm after all the other parameters since it depends on some of them (like pfc,subopt and allowlonelybasepairs)
        # Custom motif variables, custom_X is for the filepaths to the .csv files, replace_X is for if the customs should append to or replace the underlying motifs from RNA3D or Rfam
        self.motif_weighting = weight
        self.motif_fraction = fraction
        self.version = version

    @property
    def temperature(self) -> float:
        return self._temperature

    @temperature.setter
    def temperature(self, temp: float) -> None:
        if not interval_check(-273.0, 100.0, temp):
            raise ValueError("Temperature has to be between -273 and 100 Kelvin.")
        self._temperature = temp

    @property
    def energy_percent(self) -> float | None:
        return self._energy_perc

    @energy_percent.setter
    def energy_percent(self, range: float | None) -> None:
        if range is not None:
            if interval_check(float128(0), float128(100), float128(range)):
                self._energy_perc = range
            else:
                raise ValueError("Energy range cannot be below 0 or above 1.")
        elif range is None:
            self._energy_perc = range

    # Set energy range for suboptimal candidate calcualtions
    @property
    def absolute_energy(self) -> float | None:
        return self._energy

    @absolute_energy.setter
    def absolute_energy(self, e: float | None):
        if e is None:
            self._energy = None
        else:
            if e >= 0:
                self._energy = e
            else:
                raise ValueError("Energy range cannot be lower than 0.")

    # Set kvalue fpr kbest and kbacktracing
    @property
    def kvalue(self) -> int:
        return self._kvalue

    @kvalue.setter
    def kvalue(self, k: int) -> None:
        if k > 0:
            self._kvalue = k
        else:
            raise ValueError("Kvalue cannot be 0 or lower.")

    @property
    def low_probability_filter(self) -> float:
        return self._low_probability_filter

    @low_probability_filter.setter
    def low_probability_filter(self, value: float) -> None:
        if interval_check(float128(0), float128(1), float128(value)):
            self._low_probability_filter = value
        else:
            raise ValueError("Probability filter cannot be below 0 or above 1")

    # Finds path to your chosen algorithm, if it does not exist i attempts to compile the algorithm
    @property
    def algorithm_path(self) -> Path:
        path = RNAMOTIFOLD_ROOT_DIR / "bin" / self.algorithm_binary
        logger.debug(f"Set algorithm path as {path}")
        return path

    # Builds the algorithm binary name from the set parameters
    @property
    def algorithm_binary(self) -> str:
        if self.algorithm == AlgorithmType.RNAMOTIALIGN:
            if self.subopt or self.pfc:
                raise NotImplementedError(
                    "RNAmotiAlign is not available with subopt or partition function implementations"
                )
            return "RNAmotiAlign"
        match (self.allowLonelyBasepairs, self.subopt, self.pfc):
            case (2, False, False):
                alg = self.algorithm + "Motmicro"
            case (2, True, False):
                alg = self.algorithm + "_motmacro_subopt"
            case (2, False, True):
                alg = self.algorithm + "_motmacro_pfc"
            case (1, False, False) | (0, False, False):
                alg = self.algorithm
            case (1, True, False) | (0, True, False):
                alg = self.algorithm + "_subopt"
            case (1, False, True) | (0, False, True):
                alg = self.algorithm + "_pfc"
            case _:
                raise ValueError(
                    "The algorithm you specified does not exist, please revisit your arguments."
                )
        logger.debug(f"Set algorithm binary as {alg}")
        return alg

    def make_main_call(self) -> list[str]:
        """Automatic call setter. The function checks the set algorithm and builds a call string based on it."""

        runtime_dictionary: dict[str, str | int | float | Path | None] = {
            "-Q": self.motif_source,
            "-b": self.motif_orientation,
            "-t": self.temperature,
            "-v": self.version,
        }
        # Only one of these will be activate at a time
        if self.subopt:
            # Ordering here is important, the last one is always used so to keep -e overwriting -c this is necessary
            runtime_dictionary["-c"] = self.energy_percent
            runtime_dictionary["-e"] = self.absolute_energy
        elif self.pfc:
            runtime_dictionary["-F"] = self.low_probability_filter
        else:
            runtime_dictionary["-k"] = self.kvalue

        if self.algorithm == AlgorithmType.RNAMOSH:
            runtime_dictionary["-q"] = self.shape_level
        if self.algorithm == AlgorithmType.RNAMOTIALIGN:
            runtime_dictionary["-W"] = self.motif_weighting
            runtime_dictionary["-D"] = self.motif_fraction
        if (
            self.allowLonelyBasepairs == AllowLonelyBP.NO_LP
            or self.allowLonelyBasepairs == AllowLonelyBP.ALL_LP
        ):
            runtime_dictionary["-u"] = self.allowLonelyBasepairs
        arguments = [str(self.algorithm_path)]
        for k, v in runtime_dictionary.items():
            if v is not None:
                arguments.extend([k, str(v)])
        logger.debug(f"Set call as {' '.join(arguments)}")
        return arguments
