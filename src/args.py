import argparse
import configparser
import logging
from pathlib import Path
from os import cpu_count, access, W_OK
from dataclasses import dataclass


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
        return Path(pathname).resolve().exists() and is_path_creatable(pathname)
    except OSError:
        print(
            f"Given output file path is neither a file nor a dictionary that the current user can edit, defaulting to outputting to stdout"
        )
        return False


class LogCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        if not isinstance(getattr(logging, value.upper(), None), int):
            raise ValueError(f"Invalid log level: {value}")
        else:
            setattr(namespace, self.dest, value)


class WorkerCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: int, option_string):
        if value > cpu_count():
            print("Given worker number exceeds detected cpu count, setting workers to cpu_count -1")
            setattr(namespace, self.dest, cpu_count() - 1)
        else:
            setattr(namespace, self.dest, value)


class ConfigCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string=None):
        if value == "":
            print(f"Using default config: {script_parameters.default_config_path}")
            setattr(namespace, self.dest, value)

        elif Path(value).is_file():
            print(f"Using user config: {value}.")
            setattr(namespace, self.dest, Path(value))
        else:
            raise FileNotFoundError(f"Could not find specified config file {value}")


class OutputFileCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        if value == "":
            print("-o is set but no output file specified, defaulting to stdout")
            setattr(namespace, self.dest, None)
        elif is_path_exists_or_creatable(value):
            if Path(value).resolve().is_file():
                print("Output file already exists, results will be appended.")
            else:
                with open(Path(value), "a+") as file:
                    pass
            setattr(namespace, self.dest, Path(value))
        else:
            raise FileNotFoundError("Given path is not a valid path.")


@dataclass
class script_parameters:
    default_config_path = Path(__file__).resolve().parent.joinpath("data", "defaults.ini")
    user_config_path = Path(__file__).resolve().parent.joinpath("config.ini")
    RNAmotiFold_path = Path(__file__).resolve().parents[1]


def get_cmdarguments() -> argparse.Namespace:
    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(script_parameters.default_config_path))
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
        help="For interactive sessions or with single sequence as input set an ID for the output. Default is N/A.",
        dest="id",
        type=str,
        default=config.get(config.default_section, "name"),
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
        "-o",
        "--output",
        help="Set file to write results to. Results will be appended to file if it already exists! Default is no file, outputs are written to stdout",
        type=str,
        default=config.get(config.default_section, "output"),
        action=OutputFileCheck,
        nargs="?",
        const="",
        dest="output",
    )
    parser.add_argument(
        "--conf",
        help="Specify a config file path, if no path is given this defaults to RNAmotiFold/src/config.ini. Defaults are set in RNAmotiFold/src/data/defaults.ini. If --conf is set other commandline arguments will be ignored.",
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
        help="Specify which algorithm should be used, prebuild choices are: RNAmotiFold, RNAmoSh and RNAmotiCes. Set RNAmoSh shape level with -q [1-5] and RNAmotiCes mode with -p [h,b,m]. Use -s to use subopt folding. --pfc activates pfc calcualtions instead of minimum free energy. Default is RNAmotiFold",
        type=str,
        choices=[
            "RNAmotiFold",
            "RNAmoSh",
            "RNAmotiCes",
        ],
        default=config.get(config.default_section, "algorithm"),
        nargs="?",
        dest="algorithm",
    )
    pfc_or_subopt.add_argument(
        "-s",
        "--subopt",
        help="Specify if subopt folding should be used. Not useable with partition function implementations. Default is off",
        action="store_true",
        default=config.getboolean(config.default_section, "subopt"),
        dest="subopt",
    )
    parser.add_argument(
        "-Q",
        "--motif_source",
        help="Specify from which database motifs should be used, 1 = BGSU, 2 = Rfam, 3 = both, 4  = custom motifs. Default is 1",
        choices=[
            1,
            2,
            3,
            4,
        ],
        type=int,
        default=config.getint(config.default_section, "motif_source"),
        dest="motif_source",
    )
    parser.add_argument(
        "-b",
        "--orientation",
        help="Specify motif orientation: 1 = 5'-> 3',  2 = 3' -> 5' or 3 = both. Default is 1. ",
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
        help="Specify k for k-best classes get classified. Default is k = 5",
        type=int,
        default=config.getint(config.default_section, "kvalue"),
        dest="kvalue",
    )
    parser.add_argument(
        "-p",
        "--hishape",
        help="Set hishape mode, default is h",
        choices=[
            "h",
            "m",
            "b",
        ],
        type=str,
        default=config.get(config.default_section, "hishape_mode"),
        dest="hishape_mode",
    )
    parser.add_argument(
        "-q",
        "--shape_level",
        help="Set shape abstraction level. Default is 3",
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
    # Energy has to be implemented with  string since it is possibly empty and configparse can't handle it being an integer cause allow_no_value sets things to an empty string instead of None for some fucking reason.
    parser.add_argument(
        "-e",
        "--energy",
        help="Specify energy range if subopt is used. Default is None (so you can actually use the -c parameters).",
        type=str,
        default=config.get(config.default_section, "energy"),
        dest="energy",
    )
    parser.add_argument(
        "--time",
        help="Activate time logging, activating this will run predictions with unix time utility. Default is off",
        action="store_true",
        default=config.getboolean(config.default_section, "time"),
        dest="time",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        help="Scale energy parameters for folding to given temperature in Celsius. Default is 37",
        type=float,
        default=config.getfloat(config.default_section, "temperature"),
        dest="temperature",
    )
    parser.add_argument(
        "-u",
        help="Allow lonely base pairs, 1 = yes, 0 = no. Default is 0",
        type=int,
        choices=[0, 1],
        dest="basepairs",
        default=config.getint(config.default_section, "basepairs"),
    )
    parser.add_argument(
        "-c",
        help="Set energy range in %%. Is overwritten if -e is set. Default is 10.0",
        type=float,
        dest="energy_percent",
        default=config.getfloat(config.default_section, "energy_percent"),
    )
    pfc_or_subopt.add_argument(
        "--pfc",
        help="If set, calculates cumulative partition function value for each class instead of default minimum free energy predictions. Can lead to long runtimes. Default is off.",
        dest="pfc",
        action="store_true",
        default=config.getboolean(config.default_section, "pfc"),
    )
    parser.add_argument(
        "--low_prob_filter",
        help="If set, classes with a probability below 0.0001 are shown in the output when using a partition function. Default is off.",
        action="store_true",
        dest="pfc_filtering",
        default=config.getboolean(config.default_section, "pfc_filtering"),
    )

    ##############The line between script arguments and class args######
    parser.add_argument(
        "-w",
        "--workers",
        help="Specify how many predictions should be done in parallel for file input. Default is os.cpu_count()-2",
        type=int,
        action=WorkerCheck,
        default=config.getint(config.default_section, "workers"),
        dest="workers",
    )
    parser.add_argument(
        "-l",
        "--loglevel",
        help="Set log level. Default is Info. Currently nothing is getting logged to be honest, this will (hopefully) change in the future.",
        action=LogCheck,
        type=str,
        default=config[config.default_section]["loglevel"],
        dest="loglevel",
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
        "-nu",
        "--no_update",
        help="Blocks sequence updating, overwrites update_algorithms. Default is False",
        default=config.getboolean(config.default_section, "no_update"),
        action="store_true",
        dest="no_update",
    )

    parser.add_argument(
        "--update_algorithms",
        help="If set algorithms are updated to use newest motif sequence versions. Can be used when manually editing motif sequence files (adding customs for example). Gets overwritten by --no_update. Default is False.",
        action="store_true",
        dest="update_algorithms",
        default=config.getboolean(config.default_section, "update_algorithms"),
    )
    args = parser.parse_args()
    if args.config is not None:
        config.read_file(open(args.config))
        return config
    else:
        return args
