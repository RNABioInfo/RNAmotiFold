import argparse
import configparser
import logging
from pathlib import Path
from src.scirpt_parameters import script_parameters as sp
from os import cpu_count


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
            print("Given worker number exceeds detected cpu count, setting workers to cpu_count -2")
            setattr(namespace, self.dest, cpu_count() - 2)
        else:
            setattr(namespace, self.dest, value)


class FileCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        try:
            if Path(value).exists():
                print("Output file already exists, results will be appended.")
                setattr(namespace, self.dest, Path(value))
            else:
                with open(Path(value), "a+") as file:
                    pass
                setattr(namespace, self.dest, Path(value))
        except FileNotFoundError as error:
            print(
                "Could not identify valid path to specified output file, please check the output file path."
            )
            raise error


def get_cmdarguments() -> argparse.Namespace:
    config_data = get_config(sp.get_conf_path())
    # Configure parser and help message
    parser = argparse.ArgumentParser(
        prog="RNAmotiFold.py",
        description="A RNA secondary structure prediction programm with multiple functionalities for your convenience",
        epilog="GONDOR CALLS FOR AID! AND ROHAN WILL ANSWER!",
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Set input for algorithm. Running RNAmotiFold with a predefined input will not start an interactive session. Input can be a filepath, an RNA sequence or a DNA sequence. DNA sequences are silently converted to RNA.",
        type=str,
        default=None,
        const=None,
        nargs="?",
        dest="input",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Set file to write results to. Results will be appended to file if it already exists!",
        type=Path,
        default=None,
        const=None,
        action=FileCheck,
        dest="output",
        nargs="?",
    )

    # Command line arguments that control which algorithm is called with which options.
    # If you add your own partition function algorithm and want the output to have probabilities be sure to add pfc at the end of the name! This tag is used to recognize partition function algorithms by the script.
    parser.add_argument(
        "-a",
        "--algorithm",
        help="Specify which algorithm should be used, prebuild choices are: motmfepretty, motpfc, motshapeX, motshapeX_pfc, mothishape, mothishape_h_pfc, mothishape_b_pfc, mothishape_m_pfc. If you want to run a mothishape you need to specify the mode with -p, if you run mothishape with pfc enter the full name (e.g. mothishape_h_pfc). Paritition function instances should be marked with pfc at the end.",
        type=str,
        choices=[
            "motmfepretty",
            "motpfc",
            "motshapeX",
            "motshapeX_pfc",
            "mothishape",
            "mothishape_h_pfc",
            "mothishape_b_pfc",
            "mothishape_m_pfc",
        ],
        default=config_data["DEFAULTS"]["algorithm"],
        nargs="?",
        dest="algorithm",
    )
    parser.add_argument(
        "-s",
        "--subopt",
        help="Specify if subopt folding should be used. Not useable with partition function implementations. Default is off",
        action="store_true",
        default=config_data.getboolean("DEFAULTS", "subopt"),
        dest="subopt",
    )
    parser.add_argument(
        "-Q",
        "--motif_source",
        help="Specify from which database motifs should be used, 1 = BGSU, 2 = Rfam, 3 = both. Default is 3",
        choices=[
            1,
            2,
            3,
        ],
        type=int,
        default=config_data.getint("DEFAULTS", "motif_source"),
        dest="motif_source",
    )
    parser.add_argument(
        "-b",
        "--orientation",
        help="Specify motif orientation: 1 = 5'-> 3',  2 = 3' -> 5' or 3 = both. Default is 3",
        choices=[
            1,
            2,
            3,
        ],
        type=int,
        default=config_data.getint("DEFAULTS", "motif_orientation"),
        dest="motif_orientation",
    )
    parser.add_argument(
        "-k",
        "--kvalue",
        help="Specify k for k-best classes get classified. Default is k = 5",
        type=int,
        default=config_data.getint("DEFAULTS", "kvalue"),
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
        default=config_data["DEFAULTS"]["hishape_mode"],
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
        default=config_data.getint("DEFAULTS", "shape_level"),
        dest="shape_level",
    )
    parser.add_argument(
        "-e",
        "--energy",
        help="Specify energy range if subopt is used. Default is 1.0",
        type=float,
        default=config_data.getfloat("DEFAULTS", "energy"),
        dest="energy",
    )
    parser.add_argument(
        "--time",
        help="Activate time logging, activating this will run predictions with unix time utility. Default is off",
        action="store_true",
        default=config_data.getboolean("DEFAULTS", "time"),
        dest="time",
    )
    parser.add_argument(
        "-t",
        "--temperature",
        help="Scale energy parameters for folding to given temperature in Celsius. Default is 37",
        type=float,
        default=config_data.getfloat("DEFAULTS", "temperature"),
        dest="temperature",
    )
    parser.add_argument(
        "-u",
        help="Allow lonely base pairs, 1 = yes, 0 = no. Default is 0",
        type=int,
        dest="basepairs",
        default=config_data.getint("DEFAULTS", "basepairs"),
    )
    parser.add_argument(
        "-c",
        help="Set energy range in %%. Is overwritten if -e is set. Default is 10.0",
        type=float,
        dest="energy_percent",
        default=config_data.getfloat("DEFAULTS", "energy_percent"),
    )
    ##############The line between script arguments and class args######
    parser.add_argument(
        "-w",
        "--workers",
        help="Specify how many predictions should be done in parallel for file input. Default is os.cpu_count()-2",
        type=int,
        action=WorkerCheck,
        default=config_data.getint("DEFAULTS", "workers"),
        dest="workers",
    )
    parser.add_argument(
        "-l",
        "--loglevel",
        help="Set log level. Default is Info",
        action=LogCheck,
        type=str,
        default=config_data["DEFAULTS"]["loglevel"],
        dest="loglevel",
    )
    parser.add_argument(
        "-v",
        "--sep",
        help="Specify separation character for output. Default is tab for human readability.",
        type=str,
        default=config_data.get("DEFAULTS", "separator"),
        dest="separator",
    )
    # Arguments for updating motif catalogue, -r activates removal of sequences from catalogue, which ones have to be specified in Motif_collection.py update function.
    parser.add_argument(
        "-fu",
        "--force_update",
        help="Force update of sequence catalogue, gets overwritten by no_update. Default is False",
        action="store_true",
        default=config_data.getboolean("DEFAULTS", "force_update"),
        dest="force_update",
    )
    parser.add_argument(
        "-nu",
        "--no_update",
        help="Block sequence updating, overwrites force_update. Default is False",
        default=config_data.getboolean("DEFAULTS", "no_update"),
        action="store_true",
        dest="no_update",
    )
    parser.add_argument(
        "-r",
        "--remove_seq",
        help="When specified you can remove specific sequences from motifs if you do not want them to be recognized. These need to be specified in Motif_collection.py update function in the for-loop. THIS IS PERMANENT UNTIL YOU UPDATE THE SEQUENCES AGAIN (with -fu or naturally through a bgsu update). By default removes UUCAA and UACG from GNRA, GUGA from UNCG. ",
        action="store_true",
        default=config_data.getboolean("DEFAULTS", "remove_bool"),
        dest="remove_bool",
    )
    parser.add_argument(
        "--low_prob_filter",
        help="If set classes with a probability below 0.0001 are shown in the output when using a partition function. Default is off",
        action="store_true",
        dest="pfc_filtering_bool",
        default=config_data.getboolean("DEFAULTS", "pfc_filtering_bool"),
    )

    args = parser.parse_args()
    return args


def get_config(path):
    conf = configparser.ConfigParser(allow_no_value=True)
    conf.read_file(open(path))
    return conf
