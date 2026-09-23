from pathlib import Path
from dataclasses import dataclass
from typing import Literal
from argparse import Namespace
from configparser import ConfigParser
from src.RNAmotiFold.input.action_overwrites import (
    OutputFileCheck,
    WorkerCheck,
    AlgorithmMatching,
)
import logging

logger = logging.getLogger(__name__)

@dataclass
class ScriptParameters:
    """Script parameter class to hold parameters for the RNAmotiFold. Can be created from argparse.Namespace or configparser.ConfigParser. 
    All parameters are optional and have default values set in defaults.ini. Mainly here to bring together arguments from argparse and config parse for compatibility with the rest of the script.
    """

    RNAmotiFold_path = Path(__file__).resolve().parents[2]
    user_config_path = Path | None
    id: str
    input: str | None
    output: Path | None
    algorithm: Literal["RNAmoSh", "RNAmotiCes", "RNAmotiFold", "RNAmotiAlign"]
    subopt: bool
    motif_source: int
    motif_orientation: Literal[1, 2, 3]
    kvalue: int
    shape_level: int
    energy: str
    temperature: float
    basepairs: Literal[0, 1, 2]
    energy_percent: float
    pfc: bool
    low_prob_filter: float
    custom_hairpins: Path | None
    custom_internals: Path | None
    custom_bulges: Path | None
    replace_hairpins: bool
    replace_internals: bool
    replace_bulges: bool
    loglevel: str
    logfile: Path | None
    workers: int
    separator: str
    update: bool
    version: str
    fast_mode: bool
    fast_mode_merge: bool
    motif_list: str
    motif_weight: float
    motif_fraction: float
    cmake_path: Path
    gapc_path: Path
    perl_path: Path

    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(0, len(v)):
            together.append("{key}={value!r}".format(key=k[i], value=v[i]))
        return f"{classname}({', '.join(together)})"

    @property
    def process_type(self):
        if self.algorithm == "RNAmotiAlign":
            return "ali"
        else:
            return "single"

    def alg_type(self):
        if self.pfc:
            return "pfc"
        elif self.algorithm == "RNAmotiAlign":
            return "ali"
        else:
            return "mfe"

    @classmethod
    def from_argparser(cls, args: Namespace):
        if args.output:
            outpath = Path(args.output)
        else:
            outpath = None
        if args.logfile:
            logfile_path = Path(args.logfile)
        else:
            logfile_path = None
        return cls(
            id=args.id,
            input=args.input,
            output=outpath,
            algorithm=args.algorithm,
            subopt=args.subopt,
            motif_source=args.motif_source,
            motif_orientation=args.motif_orientation,
            kvalue=args.kvalue,
            shape_level=args.shape_level,
            energy=args.energy,
            temperature=args.temperature,
            basepairs=args.basepairs,
            energy_percent=args.energy_percent,
            pfc=args.pfc,
            low_prob_filter=args.low_prob_filter,
            custom_hairpins=args.custom_hairpins,
            custom_internals=args.custom_internals,
            custom_bulges=args.custom_bulges,
            replace_hairpins=args.replace_hairpins,
            replace_internals=args.replace_internals,
            replace_bulges=args.replace_bulges,
            loglevel=args.loglevel,
            logfile=logfile_path,
            workers=args.workers,
            separator=args.separator,
            update=args.update,
            version=args.version,
            fast_mode=args.fast_mode,
            fast_mode_merge=args.merge,
            motif_list=args.motif_list,
            motif_weight=args.motif_weight,
            motif_fraction=args.motif_fraction,
            cmake_path=args.cmake_path,
            gapc_path=args.gapc_path,
            perl_path=args.perl_path,
        )

    @classmethod
    def from_configparser(cls, confs: ConfigParser, section_name: str = "VARIABLES"):
        confs.set(
            section_name,
            "algorithm",
            AlgorithmMatching.algorithm_matching_function(confs.get(section_name, "algorithm")),
        )
        confs.set(
            section_name,
            "output",
            OutputFileCheck.output_file_check_function(confs.get(section_name, "output")),
        )
        confs.set(
            section_name,
            "logfile",
            OutputFileCheck.output_file_check_function(confs.get(section_name, "logfile")),
        )
        confs.set(
            section_name,
            "workers",
            str(WorkerCheck.worker_check_function(confs.get(section_name, "workers"))),
        )
        if confs.get(section_name, "output"):
            outpath = Path(confs.get(section_name, "output"))
        else:
            outpath = None
        if confs.get(section_name, "logfile"):
            logpath = Path(confs.get(section_name, "logfile"))
        else:
            logpath = None
        logger.debug(f"Read config as {confs}")
        return cls(
            id=confs.get(section_name, "id"),
            input=confs.get(section_name, "input"),
            output=outpath,
            algorithm=AlgorithmMatching.algorithm_matching_function(
                confs.get(section_name, "algorithm")
            ),
            subopt=confs.getboolean(section_name, "subopt"),
            motif_source=confs.getint(section_name, "motif_source"),
            motif_orientation=confs.getint(section_name, "motif_orientation"),  # type: ignore Idk how to "get literal" but it is in the conf and arg checks for these to only be 1,2,3
            kvalue=confs.getint(section_name, "kvalue"),
            shape_level=confs.getint(section_name, "shape_level"),
            energy=confs.get(section_name, "energy"),
            temperature=confs.getfloat(section_name, "temperature"),
            basepairs=confs.getint(section_name, "basepairs"),  # type: ignore Idk how to "get literal" but it is in the conf and arg checks for these to only be 1,2,3
            energy_percent=confs.getfloat(section_name, "energy_percent"),
            pfc=confs.getboolean(section_name, "pfc"),
            low_prob_filter=confs.getfloat(section_name, "low_prob_filter"),
            custom_hairpins=Path(confs.get(section_name, "custom_hairpins")),
            custom_internals=Path(confs.get(section_name, "custom_internals")),
            custom_bulges=Path(confs.get(section_name, "custom_bulges")),
            replace_hairpins=confs.getboolean(section_name, "replace_hairpins"),
            replace_internals=confs.getboolean(section_name, "replace_internals"),
            replace_bulges=confs.getboolean(section_name, "replace_bulges"),
            loglevel=confs.get(section_name, "loglevel"),
            logfile=logpath,
            workers=confs.getint(section_name, "workers"),
            separator=confs.get(section_name, "separator"),
            update=confs.getboolean(section_name, "update"),
            version=confs.get(section_name, "version"),
            fast_mode=confs.getboolean(section_name, "fast_mode"),
            fast_mode_merge=confs.getboolean(section_name, "merge"),
            motif_list=confs.get(section_name, "motif_list"),
            motif_weight=confs.getfloat(section_name, "motif_weight"),
            motif_fraction=confs.getfloat(section_name, "motif_fraction"),
            cmake_path=Path(confs.get("INSTALLATION", "cmake_path")),
            gapc_path=Path(confs.get("INSTALLATION", "gapc_path")),
            perl_path=Path(confs.get("INSTALLATION", "perl_path")),
        )
