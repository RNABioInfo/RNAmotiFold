import argparse
import configparser
import logging
from pathlib import Path
from os import cpu_count, access, W_OK
from dataclasses import dataclass
import sys
from typing import Any, Optional, Sequence, Literal

def is_path_creatable(pathname: str) -> bool:
    """
    `True` if the current user has sufficient permissions to create the passed
    pathname; `False` otherwise.
    """
    # Parent directory of the passed path. If empty, we substitute the current
    # working directory (CWD) instead.
    dirname = Path(pathname).resolve().parent or Path.cwd()
    return access(dirname, W_OK)


def is_path_exists_or_creatable(pathname: str) -> bool:
    """
    `True` if the passed pathname is a valid pathname for the current OS _and_
    either currently exists or is hypothetically creatable; `False` otherwise.

    This function is guaranteed to _never_ raise exceptions.
    """
    try:
        # To prevent "os" module calls from raising undesirable exceptions on
        # invalid pathnames, is_pathname_valid() is explicitly called first.
        return Path(pathname).resolve().parent.exists() and is_path_creatable(pathname)
    except OSError:
        print(
            f"Given output file path is neither a file nor a dictionary that the current user can edit, defaulting to outputting to stdout"
        )
        return False


class MotifFileCkeck(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value:Optional[str]|Sequence[Any], option_string:Optional[str] = None):
        setattr(namespace, self.dest, MotifFileCheckFunction(str(value)))


def MotifFileCheckFunction(value: Optional[str]) -> str:
    if value is not None and value != "":
        if Path(value).resolve().is_file():
            return value
        raise FileNotFoundError("Could not find specified file.")
    raise ValueError("No motif file specified")


class LogCheck(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value:Optional[str]|Sequence[Any], option_string:Optional[str] = None):
        setattr(namespace, self.dest, LogCheckFunction(str(value)))


def LogCheckFunction(value: str) -> str:
    if not isinstance(getattr(logging, value.upper(), None), int):
        raise ValueError(f"Invalid log level: {value}")
    else:
        return value.upper()


class WorkerCheck(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value:Optional[str|Sequence[Any]], option_string:Optional[str] = None):
        setattr(namespace, self.dest, WorkerCheckFunction(value)) #type:ignore 


def WorkerCheckFunction(value: Optional[str]) -> Optional[int]:
    cpus:int|None = cpu_count()
    if cpus is not None:
        if value is not None:
            if int(value) > cpus:
                print("Given worker number exceeds detected cpu count, setting workers to cpu_count - 1")
                return int(cpus - 1)
            else:
                return int(value)
        else:
            return None
    print("Could not count cpus, playing it safe and setting CPU_count to 1")
    return 1


class ConfigCheck(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value: str | Sequence[Any] | None, option_string:Optional[str]=None):
        if value == "":
            sys.stderr.write(
                f"Using default config: {script_parameters.defaults_config_path}"
            )  # Make this into a log, no need to print
            setattr(namespace, self.dest, value)
        elif Path(str(value)).resolve().is_file():
            setattr(namespace, self.dest, Path(str(value)))
        else:
            raise FileNotFoundError(f"Could not find specified config file {value}")


class OutputFileCheck(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value: str| Sequence[Any] | None, option_string:Optional[str]=None):
        setattr(namespace, self.dest, OutputFileCheckFunction(value)) #type:ignore


def OutputFileCheckFunction(value: Optional[str]):
    if value is None:
        return None
    elif is_path_exists_or_creatable(value):
        if Path(value).resolve().is_file():
            sys.stderr.write(f"Given file {value} already exists, results will be appended.\n")
        return value
    else:
        raise FileNotFoundError("Given path is not a valid path.")


class AlgorithmMatching(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value: str |Sequence[Any] | None, option_string:Optional[str] = None):
        setattr(namespace, self.dest, AlgorithmMatchingFunction(value)) #type:ignore , ignored cause of the base value typing


def AlgorithmMatchingFunction(value: str) -> Literal["RNAmoSh","RNAmotiCes","RNAmotiFold"]:
    match value.strip().lower():
            case "rnamosh":
                return "RNAmoSh"
            case "rnamotices":
                return "RNAmotiCes"
            case "rnamotifold":
                return "RNAmotiFold"
            case _:
                raise ValueError(f"Invalid algorithm specified: {value}. Valid choices are RNAmoSh, RNAmotiCes and RNAmotiFold")


def config_check(parser: configparser.ConfigParser, section_name: str = "VARIABLES"):
    parser.set(
        section_name, "algorithm", AlgorithmMatchingFunction(parser.get(section_name, "algorithm"))
    )
    parser.set(section_name, "output", OutputFileCheckFunction(parser.get(section_name, "output")))
    parser.set(section_name, "logfile", OutputFileCheckFunction(parser.get(section_name, "logfile")))
    parser.set(section_name, "workers", str(WorkerCheckFunction(parser.get(section_name, "workers"))))
    return True


@dataclass
class script_parameters:
    defaults_config_path = Path(__file__).resolve().parent.joinpath("data", "defaults.ini")
    user_config_path = Path(__file__).resolve().parent.joinpath("config.ini")
    RNAmotiFold_path = Path(__file__).resolve().parents[1]
    id: str
    input: Optional[str]
    output: Optional[Path]
    algorithm: Literal["RNAmoSh","RNAmotiCes","RNAmotiFold"]
    subopt: bool
    motif_source: int
    motif_orientation: int
    kvalue: int
    shape_level: int
    energy: str
    temperature: float
    basepairs: bool
    energy_percent: float
    pfc: bool
    low_prob_filter: float
    custom_hairpins: str
    custom_internals: str
    custom_bulges: str
    replace_hairpins: bool
    replace_internals: bool
    replace_bulges: bool
    loglevel: str
    logfile: Optional[Path]
    workers: int
    separator: str
    no_update: bool
    version: str

    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together:list[str] = []
        for i in range(0, len(v)):
            together.append("{key}={value!r}".format(key=k[i], value=v[i]))
        return f"{classname}({', '.join(together)})"

    @classmethod
    def from_argparser(cls, args: argparse.Namespace):
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
            no_update=args.no_update,
            version=args.version,
        )

    @classmethod
    def from_configparser(cls, confs: configparser.ConfigParser):
        if confs.get("VARIABLES", "output"):
            outpath = Path(confs.get("VARIABLES", "output"))
        else:
            outpath = None
        if confs.get("VARIABLES", "logfile"):
            logpath = Path(confs.get("VARIABLES", "logfile"))
        else:
            logpath = None
        return cls(
            id=confs.get("VARIABLES", "id"),
            input=confs.get("VARIABLES", "input"),
            output=outpath,
            algorithm=AlgorithmMatchingFunction(confs.get("VARIABLES", "algorithm")),
            subopt=confs.getboolean("VARIABLES", "subopt"),
            motif_source=confs.getint("VARIABLES", "motif_source"),
            motif_orientation=confs.getint("VARIABLES", "motif_orientation"),
            kvalue=confs.getint("VARIABLES", "kvalue"),
            shape_level=confs.getint("VARIABLES", "shape_level"),
            energy=confs.get("VARIABLES", "energy"),
            temperature=confs.getfloat("VARIABLES", "temperature"),
            basepairs=confs.getboolean("VARIABLES", "basepairs"),
            energy_percent=confs.getfloat("VARIABLES", "energy_percent"),
            pfc=confs.getboolean("VARIABLES", "pfc"),
            low_prob_filter=confs.getfloat("VARIABLES", "low_prob_filter"),
            custom_hairpins=confs.get("VARIABLES", "custom_hairpins"),
            custom_internals=confs.get("VARIABLES", "custom_internals"),
            custom_bulges=confs.get("VARIABLES", "custom_bulges"),
            replace_hairpins=confs.getboolean("VARIABLES", "replace_hairpins"),
            replace_internals=confs.getboolean("VARIABLES", "replace_internals"),
            replace_bulges=confs.getboolean("VARIABLES", "replace_bulges"),
            loglevel=confs.get("VARIABLES", "loglevel"),
            logfile=logpath,
            workers=confs.getint("VARIABLES", "workers"),
            separator=confs.get("VARIABLES", "separator"),
            no_update=confs.getboolean("VARIABLES", "no_update"),
            version=confs.get("VARIABLES", "version"),
        )


def get_cmdarguments() -> script_parameters:
    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(script_parameters.defaults_config_path))
    ###workaround for allow_no_value setting "option = " to an empty string (which makes sense it's just inconvenient cause it looks weird in the defaults file)
    for option in [
        x for x in config[config.default_section] if config[config.default_section][x] == ""
    ]:
        config.set(config.default_section, option, None)
    # Configure parser and help message
    parser = argparse.ArgumentParser(
        prog="RNAmotiFold.py",
        description="A RNA secondary structure prediction programm with multiple functionalities for your convenience. Starting the algorithm without an input starts an interactive session, which can be ended by inputting Exit. No interactive session will be started if you specify a RNA/DNA sequence or filepath with -i",
        epilog="GONDOR CALLS FOR AID! AND ROHAN WILL ANSWER!",
    )
    pfc_or_subopt = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "-n",
        "--name",
        help=f"For interactive sessions or with single sequence as input set an ID for the output. Default is {config.get(config.default_section, "id")}",
        dest="id",
        type=str,
        default=config.get(config.default_section, "id"),
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Set input for algorithm. Running RNAmotiFold with a predefined input will not start an interactive session. Input can be a filepath, an RNA sequence or a DNA sequence. DNA sequences are silently converted to RNA.",
        type=str,
        default=config.get(config.default_section, "input"),
        nargs="?",
        dest="input",
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Set file to write results to. Results will be appended to file if it already exists! Default is stdout.",
        type=str,
        default=config.get(config.default_section, "output"),
        action=OutputFileCheck,
        nargs="?",
        dest="output",
    )
    parser.add_argument(
        "--conf",
        help=f"Specify a config file path, if no path is given this defaults to the prewritten config file in {script_parameters.user_config_path}. Defaults are set in RNAmotiFold/src/data/defaults.ini. If --conf is set other commandline arguments will be ignored.",
        type=str,
        action=ConfigCheck,
        const=script_parameters.user_config_path,
        nargs="?",
        dest="config",
    )
    # Command line arguments that control which algorithm is called with which options.
    # If you add your own partition function algorithm and want the output to have probabilities be sure to add pfc at the end of the name! This tag is used to recognize partition function algorithms by the script.
    parser.add_argument(
        "-a",
        "--algorithm",
        help=f"Specify which algorithm should be used, prebuild choices are: RNAmotiFold, RNAmoSh and RNAmotiCes. Set RNAmoSh shape level with -q [1-5].. Use -s to use subopt folding. --pfc activates pfc calcualtions instead of minimum free energy. Default is {config.get(config.default_section, "algorithm")}",
        type=str,
        action=AlgorithmMatching,
        default=config.get(config.default_section, "algorithm"),
        nargs="?",
        dest="algorithm",
    )
    parser.add_argument(
        "-v",
        "--version",
        help=f"Specify which RNA 3D Motif sequence version you want to use. Default is the newest version. Use --no_update to disabled checking for new motif versions.",
        dest="version",
        type=str,
        default="current",
    )
    pfc_or_subopt.add_argument(
        "-s",
        "--subopt",
        help=f"Specify if subopt folding should be used. Not useable with partition function implementations. Default is {config.get(config.default_section, "subopt")}",
        action="store_true",
        default=config.getboolean(config.default_section, "subopt"),
        dest="subopt",
    )
    parser.add_argument(
        "-Q",
        "--motif_source",
        help=f"Specify from which database motifs should be used, 1 = RNA 3D Motif Atlas, 2 = Rfam, 3 = both. Default is {config.get(config.default_section, "motif_source")}",
        choices=[
            1,
            2,
            3,
        ],
        type=int,
        default=config.getint(config.default_section, "motif_source"),
        dest="motif_source",
    )
    parser.add_argument(
        "-b",
        "--orientation",
        help=f"Specify motif orientation: 1 = 5'-> 3',  2 = 3' -> 5' or 3 = both. Default is {config.get(config.default_section, "motif_orientation")}.",
        choices=[
            1,
            2,
            3,
        ],
        type=int,
        default=config.getint(config.default_section, "motif_orientation"),
        dest="motif_orientation",
    )
    parser.add_argument(
        "-k",
        "--kvalue",
        help=f"Specify k to classify only the k lowest free energy classes. Default is {config.get(config.default_section, "kvalue")}.",
        type=int,
        default=config.getint(config.default_section, "kvalue"),
        dest="kvalue",
    )
    parser.add_argument(
        "-q",
        "--shape_level",
        help=f"Set shape abstraction level. Default is {config.get(config.default_section, "shape_level")}.",
        choices=[
            1,
            2,
            3,
            4,
            5,
        ],
        type=int,
        default=config.getint(config.default_section, "shape_level"),
        dest="shape_level",
    )
    # Energy has to be implemented with  string since it is possibly empty and configparse can't handle it being an integer cause allow_no_value sets things to an empty string.
    parser.add_argument(
        "-e",
        "--energy",
        help="Specify energy range if subopt is used. Default is None (so you can actually use the -c parameters).",
        type=str,
        default=config.get(config.default_section, "energy"),
        dest="energy",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        help=f"Scale energy parameters for folding to given temperature in Celsius. Default is {config.get(config.default_section, "temperature")} °C.",
        type=float,
        default=config.getfloat(config.default_section, "temperature"),
        dest="temperature",
    )
    parser.add_argument(
        "-u",
        help=f"Allow lonely base pairs, True = yes, False = no. Default is {config.get(config.default_section, "basepairs")}",
        dest="basepairs",
        action="store_true",
        default=config.getboolean(config.default_section, "basepairs"),
    )
    parser.add_argument(
        "-c",
        help=f"Set energy range in %%. Gets overruled by -e. Default is {config.get(config.default_section, "energy_percent")}",
        type=float,
        dest="energy_percent",
        default=config.getfloat(config.default_section, "energy_percent"),
    )
    pfc_or_subopt.add_argument(
        "--pfc",
        help=f"If set, calculates cumulative partition function value for each class instead of default minimum free energy predictions. Can lead to long runtimes. Default is {config.get(config.default_section, "pfc")}.",
        dest="pfc",
        action="store_true",
        default=config.getboolean(config.default_section, "pfc"),
    )
    parser.add_argument(
        "--low_prob_filter",
        help=f"Set probability cutoff for partition function, filters out results with lower probability during calculation. Default is {config.get(config.default_section, "low_prob_filter")}.",
        type=float,
        dest="low_prob_filter",
        default=config.getfloat(config.default_section, "low_prob_filter"),
    )
    parser.add_argument(
        "-X",
        "--custom_hairpins",
        dest="custom_hairpins",
        help="Specify path to custom hairpin motif sequence csv file. File format: [sequence],[abbreviation][newline]. Check the CSV files in RNAmotiFold/src/data/motifs/ for examples.",
        action=MotifFileCkeck,
        default=config.get(config.default_section, "custom_hairpins"),
    )
    parser.add_argument(
        "-Y",
        "--custom_internals",
        dest="custom_internals",
        help="Specify path to custom internal motif sequence csv file. File format: [sequenceA]$[sequenceB],[abbreviation][newline]. Check the CSV files in RNAmotiFold/src/data/motifs/ for examples.",
        action=MotifFileCkeck,
        default=config.get(config.default_section, "custom_internals"),
    )
    parser.add_argument(
        "-Z",
        "--custom_bulges",
        dest="custom_bulges",
        help="Specify path to custom bulge motif sequence csv file. File format: [sequence],[abbreviation][newline]. Check the CSV files in RNAmotiFold/src/data/motifs/ for examples.",
        action=MotifFileCkeck,
        default=config.get(config.default_section, "custom_bulges"),
    )
    parser.add_argument(
        "-L",
        "--replace_hairpins",
        dest="replace_hairpins",
        help=f"If set, instead of appending custom hairpins to the chosen RNA 3D Motif Atlas/Rfam sequences they will fully replace them. Default is {config.get(config.default_section,"replace_hairpins")}",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_hairpins"),
    )
    parser.add_argument(
        "-E",
        "--replace_internals",
        dest="replace_internals",
        help=f"If set, instead of appending custom internals to the chosen RNA 3D Motif Atlas/Rfam sequences they will fully replace them. Default is {config.get(config.default_section,"replace_internals")}",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_internals"),
    )
    parser.add_argument(
        "-G",
        "--replace_bulges",
        dest="replace_bulges",
        help=f"If set, instead of appending custom bulges to the chosen RNA 3D Motif Atlas/Rfam sequences they will fully replace them. Default is {config.get(config.default_section,"replace_bulges")}",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_bulges"),
    )
    ##############The line between script arguments and class args######
    parser.add_argument(
        "-w",
        "--workers",
        help=f"Specify how many predictions should be done in parallel for file input. Default is {config.get(config.default_section, "workers")}",
        type=int,
        action=WorkerCheck,
        default=config.getint(config.default_section, "workers"),
        dest="workers",
    )
    parser.add_argument(
        "--loglevel",
        help=f"Set log level. Default is {config.get(config.default_section,"loglevel")}.",
        action=LogCheck,
        type=str,
        default=config.get(config.default_section, "loglevel"),
        dest="loglevel",
    )
    parser.add_argument(
        "--logfile",
        help=f"Set filepath as destination for log entries. Default is stderr stream.",
        type=str,
        action=OutputFileCheck,
        default=config.get(config.default_section, "logfile"),
        dest="logfile",
    )
    parser.add_argument(
        "--sep",
        help="Specify separation character for output. Default is tab for human readability.",
        type=str,
        default=config.get(config.default_section, "separator"),
        dest="separator",
    )
    # Arguments for updating motifs
    parser.add_argument(
        "--nu",
        "--no_update",
        help=f"Blocks checking for new RNA 3D Motif Atlas version. Default is {config.get(config.default_section, "no_update")}",
        default=config.getboolean(config.default_section, "no_update"),
        action="store_true",
        dest="no_update",
    )

    args = parser.parse_args()
    if args.config is not None:
        config.read_file(open(args.config))
        config_check(config)
        return script_parameters.from_configparser(config)
    else:
        return script_parameters.from_argparser(args)
