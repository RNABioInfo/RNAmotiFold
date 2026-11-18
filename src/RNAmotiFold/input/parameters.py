import logging
from argparse import Namespace
from configparser import ConfigParser
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import ClassVar, Self

from Bio import SeqRecord
from Bio.Align import MultipleSeqAlignment
from RNAmotiFold.input.action_overwrites import (
    AlgorithmMatching,
    OutputFileCheck,
)
from RNAmotiFold.Misc.EnumTypes import *
from RNAmotiFold.Misc.ValueObjects import *

logger = logging.getLogger(__name__)


@dataclass
class InstallParameters:
    gapcM_path: Path
    perl_path: Path
    make_path: Path
    workers: int
    force: bool

    @classmethod
    def from_namespace(cls, args: Namespace) -> Self:
        return cls(
            gapcM_path=args.gapcM_path,
            perl_path=args.perl_path,
            make_path=args.make_path,
            workers=args.workers,
            force=args.force,
        )


@dataclass
class ScriptParameters:
    """Script parameter class to hold parameters for the RNAmotiFold. Can be created from argparse.Namespace or configparser.ConfigParser.
    All parameters are optional and have default values set in defaults.ini. Mainly here to bring together arguments from argparse and config parse for compatibility with the rest of the script.
    """

    RNAmotiFold_path: ClassVar[Path] = Path(__file__).resolve().parents[2]
    id: str
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path | None
    output: Path | None
    algorithm: AlgorithmType
    subopt: bool
    motif_source: MotifSource
    motif_orientation: MotifOrientation
    kvalue: pos_int
    shape_level: ShapeLevel
    subopt_opts: SubOptEnergy
    temperature: limited_temperature
    basepairs: AllowLonelyBP
    pfc: bool
    low_prob_filter: fraction
    custom_hairpins: Path | None
    custom_internals: Path | None
    custom_bulges: Path | None
    replace_hairpins: bool
    replace_internals: bool
    replace_bulges: bool
    loglevel: str
    logfile: Path | None
    workers: pos_int
    separator: str
    version: str
    single_motif_mode: bool
    motif_list: str
    motif_weight: Decimal
    motif_fraction: fraction

    def __repr__(self) -> str:
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(len(v)):
            together.append(f"{k[i]}={v[i]!r}")
        return f"{classname}({', '.join(together)})"

    @property
    def process_type(self) -> ProcessType:
        if self.algorithm == AlgorithmType.RNAMOTIALIGN:
            return ProcessType.ALIGNMENT
        else:
            return ProcessType.SINGLE

    @property
    def alg_type(self) -> PredictionType:
        if self.algorithm == AlgorithmType.RNAMOTIALIGN:
            return PredictionType.ALI
        elif self.pfc:
            return PredictionType.PFC
        else:
            return PredictionType.MFE

    @classmethod
    def from_argparser(cls, args: Namespace) -> Self:
        if args.output:
            outpath = Path(args.output)
        else:
            outpath = None
        if args.logfile:
            logfile_path = Path(args.logfile)
        else:
            logfile_path = None
        if args.energy is not None:
            subopt_set = SubOptEnergy(pos_decimal(args.energy), SubOptMode.ABSOLUTE)
        elif args.energy_percent is not None:
            subopt_set = SubOptEnergy(
                pos_decimal(args.energy_percent), SubOptMode.PERCENT
            )
        else:
            subopt_set = SubOptEnergy(pos_decimal(10.0), SubOptMode.PERCENT)

        return cls(
            id=args.id,
            input=args.input,
            output=outpath,
            algorithm=AlgorithmType(args.algorithm),
            subopt=args.subopt,
            subopt_opts=subopt_set,
            motif_source=MotifSource(args.motif_source),
            motif_orientation=MotifOrientation(args.motif_orientation),
            kvalue=args.kvalue,
            shape_level=ShapeLevel(args.shape_level),
            temperature=args.temperature,
            basepairs=AllowLonelyBP(args.basepairs),
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
            version=args.version,
            single_motif_mode=args.single_motif_mode,
            motif_list=args.motif_list,
            motif_weight=args.motif_weight,
            motif_fraction=args.motif_fraction,
        )

    @classmethod
    def from_configparser(
        cls, confs: ConfigParser, section_name: str = "VARIABLES"
    ) -> Self:
        confs.set(
            section_name,
            "algorithm",
            AlgorithmMatching.algorithm_matching_function(
                confs.get(section_name, "algorithm")
            ),
        )
        confs.set(
            section_name,
            "output",
            OutputFileCheck.output_file_check_function(
                confs.get(section_name, "output")
            ),
        )
        confs.set(
            section_name,
            "logfile",
            OutputFileCheck.output_file_check_function(
                confs.get(section_name, "logfile")
            ),
        )
        confs.set(
            section_name,
            "workers",
            str(pos_int(confs.get(section_name, "workers"))),
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

        if confs.get(section_name, "custom_hairpins"):
            custom_hairpins = Path(confs.get(section_name, "custom_hairpins"))
        else:
            custom_hairpins = None
        if confs.get(section_name, "custom_internals"):
            custom_internals = Path(confs.get(section_name, "custom_internals"))
        else:
            custom_internals = None
        if confs.get(section_name, "custom_bulges"):
            custom_bulges = Path(confs.get(section_name, "custom_bulges"))
        else:
            custom_bulges = None
        if confs.get(section_name, "energy"):
            subopt_set = SubOptEnergy(
                pos_decimal(confs.get(section_name, "energy")),
                SubOptMode.ABSOLUTE,
            )
        elif confs.get(section_name, "energy_percent"):
            subopt_set = SubOptEnergy(
                pos_decimal(confs.get(section_name, "energy_percent")),
                SubOptMode.PERCENT,
            )
        else:
            subopt_set = SubOptEnergy(pos_decimal(10.0), SubOptMode.PERCENT)

        return cls(
            id=confs.get(section_name, "id"),
            input=confs.get(section_name, "input"),
            output=outpath,
            algorithm=AlgorithmMatching.algorithm_matching_function(
                confs.get(section_name, "algorithm")
            ),
            subopt=confs.getboolean(section_name, "subopt"),
            motif_source=MotifSource(confs.getint(section_name, "motif_source")),
            motif_orientation=confs.getint(section_name, "motif_orientation"),  # type: ignore Idk how to "get literal" but it is in the conf and arg checks for these to only be 1,2,3
            kvalue=pos_int(confs.getint(section_name, "kvalue")),
            shape_level=ShapeLevel(confs.getint(section_name, "shape_level")),
            subopt_opts=subopt_set,
            temperature=limited_temperature(confs.get(section_name, "temperature")),
            basepairs=confs.getint(section_name, "basepairs"),  # type: ignore Idk how to "get literal" but it is in the conf and arg checks for these to only be 1,2,3
            pfc=confs.getboolean(section_name, "pfc"),
            low_prob_filter=fraction(confs.getfloat(section_name, "low_prob_filter")),
            custom_hairpins=custom_hairpins,
            custom_internals=custom_internals,
            custom_bulges=custom_bulges,
            replace_hairpins=confs.getboolean(section_name, "replace_hairpins"),
            replace_internals=confs.getboolean(section_name, "replace_internals"),
            replace_bulges=confs.getboolean(section_name, "replace_bulges"),
            loglevel=confs.get(section_name, "loglevel"),
            logfile=logpath,
            workers=pos_int(confs.getint(section_name, "workers")),
            separator=confs.get(section_name, "separator"),
            version=confs.get(section_name, "version"),
            single_motif_mode=confs.getboolean(section_name, "single_motif_mode"),
            motif_list=confs.get(section_name, "motifs"),
            motif_weight=Decimal(confs.getfloat(section_name, "motif_weight")),
            motif_fraction=fraction(confs.getfloat(section_name, "motif_fraction")),
        )
