import argparse
import configparser
import logging
from pathlib import Path
from src.input.parameters import ScriptParameters
from src.input.action_overwrites import (
    MotifFileCheck,
    LogCheck,
    FloatCheck,
    WorkerCheck,
    OutputFileCheck,
    ConfigCheck,
    MotifListCheck,
    AlgorithmMatching,
)

loggers = logging.getLogger("RNAmotiFold.args")
_defaults_config_path = Path(__file__).resolve().parents[1].joinpath("defaults", "defaults.ini")


def get_cmdarguments() -> tuple[ScriptParameters, list[str]]:
    """Sets up argument parser using defaults/defaults.ini for default values. Checks set arguments using action_overwrites and will raise Errors if something is not right. Returns a tuple of ScriptParameters and a list of unknown arguments"""
    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(_defaults_config_path))
    ###workaround for allow_no_value setting "option = " to an empty string (which makes sense it's just inconvenient cause it looks weird in the defaults file)
    for option in [
        x for x in config[config.default_section] if config[config.default_section][x] == ""
    ]:
        config.set(config.default_section, option, None)
    # Configure parser and help message
    parser = argparse.ArgumentParser(
        prog="RNAmotiFold.py",
        description="A RNA secondary structure prediction programm with multiple functionalities for your convenience. Starting the algorithm without an input starts an interactive session, which can be ended by inputting Exit. No interactive session will be started if you specify a RNA/DNA sequence or filepath with -i. Defaults are set in RNAmotiFold/src/data/defaults.ini and count for both config files and commandline arguments.",
        epilog="GONDOR CALLS FOR AID! AND ROHAN WILL ANSWER!",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    pfc_or_subopt = parser.add_mutually_exclusive_group()
    parser.add_argument(
        "-n",
        "--name",
        help=f"For interactive sessions or with single sequence as input set an ID for the output.",
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
        help="Set file to write results to. Results will be appended to file if it already exists! If set to None results will be printed to stdout.",
        type=str,
        default=config.get(config.default_section, "output"),
        action=OutputFileCheck,
        nargs="?",
        dest="output",
    )
    parser.add_argument(
        "--conf",
        help=f"Specify a config file path, if no path is given this defaults to the prewritten config file at {ScriptParameters.user_config_path}.  If --conf is set other commandline arguments will be ignored.",
        type=str,
        action=ConfigCheck,
        const=ScriptParameters.user_config_path,
        nargs="?",
        dest="config",
    )
    parser.add_argument(
        "-f",
        "--single-motif",
        default=config.getboolean(config.default_section, "fast_mode"),
        action="store_true",
        dest="fast_mode",
        help=f"Enables single-motif mode rediction mode, instead of all motifs being predicted at once they are each predicted separately and merged afterwards.",
    )
    parser.add_argument(
        "-m",
        "--fast_mode_merge",
        action="store_true",
        default=config.getboolean(config.default_section, "fast_mode_merge"),
        dest="merge",
        help="Enables merging of structures during fast mode, allowing for combined motif outputs even with fast mode. Currently only implemented with RNAmotiFold!",
    )
    # Command line arguments that control which algorithm is called with which options.
    # If you add your own partition function algorithm and want the output to have probabilities be sure to add pfc at the end of the name! This tag is used to recognize partition function algorithms by the script.
    parser.add_argument(
        "-a",
        "--algorithm",
        help=f"Specify which algorithm should be used, prebuild choices are: RNAmotiFold, RNAmoSh and RNAmotiCes. Set RNAmoSh shape level with -q [1-5].. Use -s to use subopt folding. --pfc activates pfc calcualtions instead of minimum free energy.",
        type=str,
        action=AlgorithmMatching,
        default=config.get(config.default_section, "algorithm"),
        nargs="?",
        dest="algorithm",
    )
    parser.add_argument(
        "-v",
        "--version",
        help=f"Specify which RNA 3D Motif sequence version you want to use. Use --no_update to disabled checking for new motif versions.",
        dest="version",
        type=str,
        default="current",
    )
    pfc_or_subopt.add_argument(
        "--s",
        "--subopt",
        help=f"Specify if subopt folding should be used. Not useable with partition function implementations.",
        action="store_true",
        default=config.getboolean(config.default_section, "subopt"),
        dest="subopt",
    )
    parser.add_argument(
        "-Q",
        "--motif_source",
        help=f"Specify from which database motifs should be used, 1 = RNA 3D Motif Atlas, 2 = Rfam, 3 = both.",
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
        help=f"Specify motif orientation: 1 = 5'-> 3',  2 = 3' -> 5' or 3 = both.",
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
        help=f"Specify k to classify only the k lowest free energy classes.",
        type=int,
        default=config.getint(config.default_section, "kvalue"),
        dest="kvalue",
    )
    parser.add_argument(
        "-q",
        "--shape_level",
        help=f"Set shape abstraction level.",
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
        help="Specify energy range if subopt is used. Defaults to None so you can actually use the -c parameters.",
        type=str,
        default=config.get(config.default_section, "energy"),
        dest="energy",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        help=f"Scale energy parameters for folding to given temperature in Celsius.",
        type=float,
        default=config.getfloat(config.default_section, "temperature"),
        dest="temperature",
    )
    parser.add_argument(
        "-u",
        "--allow_lonely_basepairs",
        help=f"Allow lonely base pairs can only be set to 0 (no lonely base pairs), 1 (allow all lonely base pairs),2 (allow lonely base pairs around motifs only)",
        dest="basepairs",
        choices=[0, 1, 2],
        type=int,
        default=config.getint(config.default_section, "basepairs"),
    )
    parser.add_argument(
        "-c",
        help=f"Set energy range in %%. Gets overruled by -e.",
        type=float,
        dest="energy_percent",
        default=config.getfloat(config.default_section, "energy_percent"),
    )
    pfc_or_subopt.add_argument(
        "--pfc",
        help=f"If set, calculates cumulative partition function value for each class instead of default minimum free energy predictions.",
        dest="pfc",
        action="store_true",
        default=config.getboolean(config.default_section, "pfc"),
    )
    parser.add_argument(
        "--low_prob_filter",
        help=f"Set probability cutoff for partition function, filters out results with lower probability during calculation.",
        type=float,
        dest="low_prob_filter",
        default=config.getfloat(config.default_section, "low_prob_filter"),
    )
    parser.add_argument(
        "-X",
        "--custom_hairpins",
        dest="custom_hairpins",
        help="Specify path to custom hairpin motif sequence csv file. File format: [sequence],[abbreviation][newline]. Check the CSV files in RNAmotiFold/src/data/motifs/ for examples.",
        action=MotifFileCheck,
        default=config.get(config.default_section, "custom_hairpins"),
    )
    parser.add_argument(
        "-Y",
        "--custom_internals",
        dest="custom_internals",
        help="Specify path to custom internal motif sequence csv file. File format: [sequenceA]$[sequenceB],[abbreviation][newline]. Check the CSV files in RNAmotiFold/src/data/motifs/ for examples.",
        action=MotifFileCheck,
        default=config.get(config.default_section, "custom_internals"),
    )
    parser.add_argument(
        "-Z",
        "--custom_bulges",
        dest="custom_bulges",
        help="Specify path to custom bulge motif sequence csv file. File format: [sequence],[abbreviation][newline]. Check the CSV files in RNAmotiFold/src/data/motifs/ for examples.",
        action=MotifFileCheck,
        default=config.get(config.default_section, "custom_bulges"),
    )
    parser.add_argument(
        "-L",
        "--replace_hairpins",
        dest="replace_hairpins",
        help=f"If set, instead of appending custom hairpins to the chosen RNA 3D Motif Atlas/Rfam sequences they will fully replace them.",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_hairpins"),
    )
    parser.add_argument(
        "-E",
        "--replace_internals",
        dest="replace_internals",
        help=f"If set, instead of appending custom internals to the chosen RNA 3D Motif Atlas/Rfam sequences they will fully replace them.",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_internals"),
    )
    parser.add_argument(
        "-G",
        "--replace_bulges",
        dest="replace_bulges",
        help=f"If set, instead of appending custom bulges to the chosen RNA 3D Motif Atlas/Rfam sequences they will fully replace them.",
        action="store_true",
        default=config.getboolean(config.default_section, "replace_bulges"),
    )
    ##############The line between script arguments and class args######
    parser.add_argument(
        "-w",
        "--workers",
        help=f"Specify how many predictions should be done in parallel for file input.",
        type=int,
        action=WorkerCheck,
        default=config.getint(config.default_section, "workers"),
        dest="workers",
    )
    parser.add_argument(
        "--loglevel",
        help=f"Set log level.",
        action=LogCheck,
        type=str,
        default=config.get(config.default_section, "loglevel"),
        dest="loglevel",
    )
    parser.add_argument(
        "--logfile",
        help=f"Set filepath as destination for log entries. If set to None error messages are printed to stderr.",
        type=str,
        action=OutputFileCheck,
        default=config.get(config.default_section, "logfile"),
        dest="logfile",
    )
    parser.add_argument(
        "--sep",
        help="Specify separation character for output.",
        type=str,
        default=config.get(config.default_section, "separator"),
        dest="separator",
    )
    # Arguments for updating motifs
    parser.add_argument(
        "--nu",
        "--no_update",
        help=f"Blocks checking for new RNA 3D Motif Atlas version (saves quite some time on startup cause the server is slow).",
        default=config.getboolean(config.default_section, "no_update"),
        action="store_true",
        dest="no_update",
    )
    parser.add_argument(
        "--motifs",
        help=f"Specify which motifs should be recognized during prediction. Works with custom motifs and all modes.",
        default=config.get(config.default_section, "motif_list"),
        type=str,
        action=MotifListCheck,
        dest="motif_list",
    )
    parser.add_argument(
        "--weight",
        help=f"Specify weighting of motifs during alignment folding, multiplies motif score by this value. Default is 1.0",
        default=config.getfloat(config.default_section, "motif_weight"),
        type=float,
        dest="motif_weight",
    )
    parser.add_argument(
        "--fraction",
        help=f"Specify in how many rows a motif has to be recognized in the same loop to be accepted. Default is 0.7",
        default=config.getfloat(config.default_section, "motif_fraction"),
        type=float,
        action=FloatCheck,
        dest="motif_fraction",
    )

    args = parser.parse_known_args()
    # Some lazily done arg checks to avoid specific arg combinations that dont work or arent implemented, clean this up at some point!
    if args[0].algorithm != "RNAmotiFold" and args[0].merge:
        raise parser.error("Fast mode merging is only implemented for RNAmotiFold, sorry!")
    if args[0].config is not None:
        config.read_file(open(args[0].config))
        return (ScriptParameters.from_configparser(config), args[1])
    else:
        return (ScriptParameters.from_argparser(args[0]), args[1])
