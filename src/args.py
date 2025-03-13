import argparse
import configparser
import logging
from pathlib import Path
from os import cpu_count, access, W_OK
from dataclasses import dataclass
import sys


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
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value, option_string):
        setattr(namespace, self.dest, MotifFileCheckFunction(value))


def MotifFileCheckFunction(value: str) -> str:
    if value is not None and value != "":
        if Path(value).resolve().is_file():
            return value
        else:
            raise FileNotFoundError("Could not find specified file.")


class LogCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        setattr(namespace, self.dest, LogCheckFunction(value))


def LogCheckFunction(value: str) -> str:
    if not isinstance(getattr(logging, value.upper(), None), int):
        raise ValueError(f"Invalid log level: {value}")
    else:
        return value.upper()


class WorkerCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: int, option_string):
        setattr(namespace, self.dest, WorkerCheckFunction(value))


def WorkerCheckFunction(value: int | str) -> int | str:
    if int(value) > cpu_count():
        print("Given worker number exceeds detected cpu count, setting workers to cpu_count - 1")
        if isinstance(value, str):
            return str(cpu_count() - 1)
        else:
            return cpu_count() - 1
    else:
        return value


class ConfigCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string=None):
        if value == "":
            sys.stderr.write(
                f"Using default config: {script_parameters.defaults_config_path}"
            )  # Make this into a log, no need to print
            setattr(namespace, self.dest, value)

        elif Path(value).resolve().is_file():
            setattr(namespace, self.dest, Path(value))
        else:
            raise FileNotFoundError(f"Could not find specified config file {value}")


class OutputFileCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        setattr(namespace, self.dest, OutputFileCheckFunction(value, self.dest))


def OutputFileCheckFunction(value: str, dest: str):
    if value is None:
        return None
    elif is_path_exists_or_creatable(value):
        if Path(value).resolve().is_file():
            sys.stderr.write(f"Given {dest} file already exists, results will be appended.\n")
        return Path(value)
    else:
        raise FileNotFoundError("Given path is not a valid path.")


class AlgorithmMatching(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        setattr(namespace, self.dest, AlgorithmMatchingFunction(value))


def AlgorithmMatchingFunction(value: str):
    match value.strip().lower():
        case "rnamosh":
            return "RNAmoSh"
        case "rnamotices":
            return "RNAmotiCes"
        case "rnamotifold":
            return "RNAmotiFold"
        case _:
            return value


def config_check(parser: configparser.ConfigParser, section_name: str = "VARIABLES"):
    parser.set(
        section_name, "algorithm", AlgorithmMatchingFunction(parser.get(section_name, "algorithm"))
    )
    parser.set(section_name, "output", OutputFileCheckFunction(parser.get(section_name, "output")))
    parser.set(section_name, "workers", WorkerCheckFunction(parser.get(section_name, "workers")))
    return True


@dataclass
class script_parameters:
    defaults_config_path = Path(__file__).resolve().parent.joinpath("data", "defaults.ini")
    user_config_path = Path(__file__).resolve().parent.joinpath("config.ini")
    RNAmotiFold_path = Path(__file__).resolve().parents[1]


def get_cmdarguments() -> argparse.Namespace:
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
        help=f"For interactive sessions or with single sequence as input set an ID for the output. Default is {config.get(config.default_section, "ID")}",
        dest="id",
        type=str,
        default=config.get(config.default_section, "ID"),
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
        help=f"Specify from which database motifs should be used, 1 = BGSU, 2 = Rfam, 3 = both. Default is {config.get(config.default_section, "motif_source")}",
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
        help=f"If set, classes with a probability below 0.0001 are shown in the output when using a partition function. Default is {config.get(config.default_section, "pfc_filtering")}.",
        action="store_true",
        dest="pfc_filtering",
        default=config.getboolean(config.default_section, "pfc_filtering"),
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
        help=f"If set, instead of appending custom hairpins to the chosen RNA3D Atlas/Rfam sequences they will fully replace them. Default is {config.get(config.default_section,"replace_hairpins")}",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_hairpins"),
    )
    parser.add_argument(
        "-E",
        "--replace_internals",
        dest="replace_internals",
        help=f"If set, instead of appending custom internals to the chosen RNA3D Atlas/Rfam sequences they will fully replace them. Default is {config.get(config.default_section,"replace_internals")}",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_internals"),
    )
    parser.add_argument(
        "-G",
        "--replace_bulges",
        dest="replace_bulges",
        help=f"If set, instead of appending custom bulges to the chosen RNA3D Atlas/Rfam sequences they will fully replace them. Default is {config.get(config.default_section,"replace_bulges")}",
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
        "-v",
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
    parser.add_argument(
        "--ua",
        "--update_algorithms",
        help=f"If set, loads current motif sequences from /RNAmotiFold/src/data/motifs/ into algorithms and recompiles them. This does not fetch the current motif sequences, re-run setup.py to fetch latest motif sequences and load them into the algorithms. Default is {config.get(config.default_section, "update_algorithms")}.",
        action="store_true",
        dest="update_algorithms",
        default=config.getboolean(config.default_section, "update_algorithms"),
    )
    args = parser.parse_args()
    if args.config is not None:
        config.read_file(open(args.config))
        config_check(config)
        # Create a log entry to record inputs
        return config
    else:
        return args
