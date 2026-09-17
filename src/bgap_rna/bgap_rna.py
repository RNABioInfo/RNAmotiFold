from pathlib import Path
from typing import Literal
import logging
from src.input.parameters import ScriptParameters

logger = logging.getLogger("bgap_rna")

# A Python class for making Bellman's GAP more convenient to use
# Just create a class instances, feed it with the call arguments you need
# and it'll create a call from the arguments classified as runtime arguments.
# These are: motif_source, motif_orientation, kvalue, hishape_mode, shape_level, energy


class bgap_rna:
    """Main class for running RNAmotiFold algortihms through python, just hand your arguments to this class (everything else will be defaults set in the class)
    and use [your class obj].auto_run([input]) to run predictions.). Most functionality is handled by the handlers (input reading, which motifs to use etc.).
    The bgap obj itself is mostly here to bring everything together. Every bgap_obj creates it's own handlers
    """

    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(0, len(v)):
            together.append("{key}={value!r}".format(key=k[i], value=v[i]))
        return f"{classname}({', '.join(together)})"

    def __str__(self) -> str:
        return self.call

    #    @staticmethod
    #    def postprocessing_mfe(merged_output: src.results.base_result.algorithm_output) -> src.results.base_result.algorithm_output:
    #        """
    #        Postprocessing function for merging outputs of the seperated motif predictions
    #        """
    #        mfe_dict: dict[float, list[base_result.result_mfe]] = {}
    #        for res in merged_output.results:
    #            if (
    #                isinstance(res, base_result.result_mfe) and res.classifier
    #            ):  # this is a little unnecessary but it gets rid of warnings, the res classifier filter removes the "no motif" structure
    #                if (
    #                    res.free_energy not in mfe_dict.keys()
    #                ):  # -> It makes no sense to have it in the merging process since if it can fit a motif it will be the mfe for that motif anyways
    #                    mfe_dict[res.free_energy] = [res]
    #                else:
    #                    mfe_dict[res.free_energy].append(res)
    #        for key in mfe_dict.keys():
    #            if len(mfe_dict[key]) > 1:
    #                merge_candidates = base_result.result_mfe.get_compatible_structures(mfe_dict[key])
    #                for compatible_structures in merge_candidates:
    #                    new_result = base_result.result_mfe.merge_structures(
    #                        [mfe_dict[key][i] for i in compatible_structures]
    #                    )
    #                    if new_result is not None:
    #                        merged_output.results.append(new_result)
    #            else:
    #                continue
    #        merged_output.results.sort(key=lambda x: x.free_energy)  # type: ignore
    #        return merged_output
    #
    #    @staticmethod
    #    def postprocessing_pfc(
    #        merged_output: list[base_result.algorithm_output],
    #    ) -> list[base_result.algorithm_output]:
    #        returnlist: list[base_result.algorithm_output] = []
    #        checklist: list[str] = []
    #        for output in merged_output:
    #            if str(output) not in checklist and len(output.results) > 1:
    #                checklist.append(str(output))
    #                returnlist.append(output)
    #        return returnlist

    @classmethod
    def from_script_parameters(cls, params: ScriptParameters):
        return cls(
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
            session_id=params.id,
            fast_mode=params.fast_mode,
            motif_string=params.motif_list,
            weight=params.motif_weight,
            fraction=params.motif_fraction,
        )

    def __init__(
        self,
        alg: Literal["RNAmotiFold", "RNAmoSh", "RNAmotiCes", "RNAmotiAlign"],
        motif_source: int = 1,
        motif_orientation: Literal[1, 2, 3] = 1,
        kvalue: int = 5,
        shape_level: int = 3,
        energy: str | None = None,
        temperature: float = 37.0,
        energy_percent: float = 5.0,
        allowLonelyBasepairs: Literal[0, 1, 2] = 0,
        subopt: bool = False,
        pfc: bool = False,
        low_prob_filter: float = 0.000001,
        session_id: str = "N/A",
        fast_mode: bool = False,
        motif_string: str = "",
        weight: float = 1.0,
        fraction: float = 0.7,
    ):
        self.id = session_id
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
        self.fast_mode = fast_mode
        self.motif_string = motif_string
        self.motif_weighting = weight
        self.motif_fraction = fraction

    # Slightly controversial addition, if custom_call is set it permanently overwrites the default call and is even returned whenever
    # the standard self.call is asked for. This avoids duplicating and overcomplicating code down the line. Deleting this will return normal calls

    @property
    def process_type(self) -> Literal["single", "ali"]:
        if self.algorithm in ["RNAmoSh", "RNAmotiCes", "RNAmotiFold"]:
            return "single"
        elif self.algorithm in ["RNAmotiAlign"]:
            return "ali"
        raise ValueError(
            "Could not identify process type as single folding or alignment folding, check your set algorithm"
        )

    @property
    def custom_call(self):
        return self._custom_call

    @custom_call.setter
    def custom_call(self, call: str):
        self._custom_call = call

    @custom_call.deleter
    def custom_call(self):
        del self._custom_call

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, temp: float):
        if not -273 < temp < 100:
            logger.info("Temperature outside realistic range, beware results may be inaccurate")
        self._temperature = temp

    @property
    def energy_percent(self):
        return self._energy_perc

    @energy_percent.setter
    def energy_percent(self, range: float | None):
        if range is not None:
            if range > 0:
                self._energy_perc = range
            else:
                raise ValueError("Energy range cannot be below 0.")
        elif range is None:
            self._energy_perc = range

    @property
    def allowLonelyBasepairs(self):
        return self._allowLonelyBasepairs

    @allowLonelyBasepairs.setter
    def allowLonelyBasepairs(self, val: Literal[0, 1, 2]):
        if val in [0, 1, 2]:
            self._allowLonelyBasepairs = val
        else:
            raise ValueError(
                "Allow lonely base pairs can only be set to 0 (no lonely base pairs), 1 (allow all lonely base pairs),2 (allow lonely base pairs around motifs only)"
            )

    # Choose motif source, 1 = RNA 3D Motif Atlas, 2 = RMFam, 3 = Both
    @property
    def motif_source(self):
        return self._motif_source

    @motif_source.setter
    def motif_source(self, Q: int):
        if Q in [1, 2, 3]:
            self._motif_source = Q
        else:
            raise ValueError(
                "Motif source can only be 1 = RNA 3D Motif Atlas , 2 = RMfam , 3 = Both"
            )

    # Choose motif orientation, 1 = 5'->3' only, 2 = 3'->5' only, 3= both
    @property
    def motif_orientation(self) -> Literal[1, 2, 3]:
        return self._motif_orientation

    @motif_orientation.setter
    def motif_orientation(self, b: Literal[1, 2, 3]):
        if b in [1, 2, 3]:
            self._motif_orientation: Literal[1] | Literal[2] | Literal[3] = b
        else:
            raise ValueError("Motif direction can only be 1 = 5'->3' , 2 = 3'->5' , 3 = Both.")

    # Set shape abstraction level, viable inputs are 1-5
    @property
    def shape_level(self):
        return self._shape_level

    @shape_level.setter
    def shape_level(self, q: int):
        if q in [1, 2, 3, 4, 5]:
            self._shape_level = q
        else:
            raise ValueError(
                "Shape level can only be set to levels 5 (most abstract) - 1 (least abstract)."
            )

    # Set energy range for suboptimal candidate calcualtions
    @property
    def absolute_energy(self):
        return self._energy

    @absolute_energy.setter
    def absolute_energy(self, e: str | None):
        if e == "":
            e = None
        if e is not None:
            if float(e) >= 0:
                self._energy = e
            else:
                raise ValueError("Energy range cannot be lower than 0.")
        else:
            self._energy = e

    # Set kvalue fpr kbest and kbacktracing
    @property
    def kvalue(self):
        return self._kvalue

    @kvalue.setter
    def kvalue(self, k: int):
        if k > 0:
            self._kvalue = k
        else:
            raise ValueError("Kvalue cannot be 0 or lower.")

    @property
    def low_probability_filter(self):
        return self._low_probability_filter

    @low_probability_filter.setter
    def low_probability_filter(self, value: float):
        if 0 <= value < 1:
            self._low_probability_filter = value
        else:
            raise ValueError("Probability filter cannot be below 0 or above 1")

    # Finds path to your chosen algorithm, if it does not exist i attempts to compile the algorithm
    @property
    def algorithm_path(self):
        return str(
            Path(__file__)
            .resolve()
            .parents[2]
            .joinpath("Build", "bin")
            .joinpath(self.algorithm_binary)
        )

    # Builds the algorithm binary name from the set parameters
    @property
    def algorithm_binary(self) -> str:
        if self.algorithm == "RNAmotiAlign":
            if self.subopt or self.pfc:
                raise NotImplementedError(
                    "RNAmotiAlign is not available with subopt or partition function implementations"
                )
            return "RNAmotiAlign"
        match (self.allowLonelyBasepairs, self.subopt, self.pfc):
            case (2, False, False):
                return self.algorithm + "Motmicro"
            case (2, True, False):
                return self.algorithm + "_motmacro_subopt"
            case (2, False, True):
                return self.algorithm + "_motmacro_pfc"
            case (1, False, False) | (0, False, False):
                return self.algorithm
            case (1, True, False) | (0, True, False):
                return self.algorithm + "_subopt"
            case (1, False, True) | (0, False, True):
                return self.algorithm + "_pfc"
            case _:
                raise ValueError(
                    "The algorithm you specified does not exist, please revisit your arguments."
                )

    @property
    def custom_bulges(self) -> Path | None:
        return self._custom_bulges

    @custom_bulges.setter
    def custom_bulges(self, path: str | Path | None) -> None:
        if path is not None:
            if Path(path).is_file():
                self._custom_bulges = Path(path)
            else:
                raise FileNotFoundError(
                    "Unabled to find specified custom bulges file, please check the file path"
                )
        else:
            self._custom_bulges = None

    @property
    def call(self):
        """Automatic call setter, if a custom_call is set this will always return the custom_call. The function checks the set algorithm and builds a call string based on it."""
        if hasattr(self, "custom_call"):
            return self.custom_call
        runtime_dictionary: dict[str, str | int | float | None | Path] = {
            "-Q": self.motif_source,
            "-b": self.motif_orientation,
            "-t": self.temperature,
        }
        if self.subopt:
            # Ordering here is important, the last one is always used so to keep -e overwriting -c this is necessary
            runtime_dictionary["-c"] = self.energy_percent
            runtime_dictionary["-e"] = self.absolute_energy
        elif self.pfc:
            runtime_dictionary["-F"] = self.low_probability_filter
        else:
            runtime_dictionary["-k"] = self.kvalue

        if self.algorithm == "RNAmoSh":
            runtime_dictionary["-q"] = self.shape_level
        if self.algorithm == "RNAmotiAlign":
            runtime_dictionary["-W"] = self.motif_weighting
            runtime_dictionary["-D"] = self.motif_fraction
        if self.allowLonelyBasepairs in [0, 1]:
            runtime_dictionary["-u"] = self.allowLonelyBasepairs
        arguments = [
            "{key} {value}".format(
                key=x,
                value=y,
            )
            for x, y in runtime_dictionary.items()
            if y is not None
        ]

        seq_free_call = " ".join(
            [self.algorithm_path, " ".join(arguments), ""]  # "/usr/bin/time"
        )  # Creates call string without a sequence
        return seq_free_call

    @call.setter
    def call(self):
        raise ValueError("Please use the custom_call property so set a custom call.")

    @property
    def motif_string(self) -> str:
        return self._motif_string

    @motif_string.setter
    def motif_string(self, motif_str: str | None) -> None:
        if motif_str is None:
            self._motif_string = ""
        else:
            self._motif_string = motif_str

    # Calibrate results.algorithm objects based on the current status of the bgap_rna class instance
#    def _calibrate_result_objects(self, sep: str = "\t"):
#        """Calibrate result objects to current configuration (separator, algorithm and pfc + probability filtering)"""
#        base_result.result.separator = sep
#        if self.pfc:
#            base_result.algorithm_output.Status = "pfc"
#        elif self.algorithm == "RNAmotiAlign":
#            base_result.algorithm_output.Status = "alignment"
#        else:
#            base_result.algorithm_output.Status = "mfe"
#
#    # Calibrate result objects and run either a single process if the input is a SeqRecord Object or run multiple predictions if input was a file or list of SeqRecord objects
#    def auto_run(
#        self,
#        user_input: (
#            SeqRecord
#            | list[SeqRecord]
#            | FastaIO.FastaIterator
#            | QualityIO.FastqPhredIterator
#            | Generator[SeqRecord, None, None]
#        ),
#        version: str,
#        o_file: Optional[Path | str] = None,
#        pool_workers: int = multiprocessing.cpu_count(),
#        output_csv_separator: str = "\t",
#        merge: bool = False,
#        name: str = "N/A",
#    ) -> list[base_result.algorithm_output | base_result.error]:
#        """Checks type of self.input and runs a Single Process in case of a SeqRecord or a MultiProcess in case of an Iterable as input."""
#
#        if output_csv_separator == r"\t":
#            output_csv_separator = output_csv_separator.replace(r"\t", "\t")
#        self._calibrate_result_objects(output_csv_separator)
#        motif_files = self._calibrate_self(version=version)
#        if self.algorithm == "RNAmotiAlign" and (
#            isinstance(user_input, FastaIO.FastaIterator) or isinstance(user_input, list)
#        ):
#            output = self._run_alignment_folding(user_input, o_file, name)
#        elif self.fast_mode:
#            output = self.run_separate_processes(
#                user_input, motif_files, o_file, pool_workers, merge
#            )
#        elif isinstance(user_input, SeqRecord):
#            output = self._run_single_process(user_input, o_file)
#        else:
#            output = self._run_multi_process(user_input, o_file, workers=pool_workers)
#        return output
