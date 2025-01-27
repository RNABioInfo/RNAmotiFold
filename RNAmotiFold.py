import src.process as proc
import src.args as args
from pathlib import Path
from argparse import Namespace
from src.results import algorithm_output
import setup
import configparser
import setup
import json
import logging


def make_new_logger(
    lvl: str | int, name: str, format: str = "", propagate: bool = False
) -> logging.Logger:
    """Convenient logger creation function,"""
    logger = logging.getLogger(name)
    logger.setLevel(lvl.upper())
    handler = logging.StreamHandler()  # logs to stderr
    if format:
        formatter = logging.Formatter(fmt=format)
    else:
        formatter = logging.Formatter(fmt="%(asctime)s:%(levelname)s: %(message)s")
    handler.setFormatter(formatter)
    # Permanent propagate False since I'm not utilizing any subloggers and sublcasses.
    logger.propagate = propagate
    logger.addHandler(handler)
    return logger


def _version_checks(force_alg_update: bool):
    defaults_config = configparser.ConfigParser(allow_no_value=True)
    defaults_config.read_file(open(args.script_parameters.default_config_path))
    try:
        update_needed = proc.bgap_rna.check_motif_versions(defaults_config["VERSIONS"]["motifs"])
        if update_needed:
            print("There is a new set of RNA 3D Motif sequences available, update ? [y/n]")
            answer = input()
            if answer in ["y", "yes", "Y", "Yes", "YEs", "YES"]:
                print("Updating...")
                setup.update_sequences_algorithms()  # Fetch newest sequences and update algorithms
                defaults_config.set("VERSIONS", "motifs", proc.bgap_rna.get_current_motif_version())
                with open(args.script_parameters.default_config_path, "w") as cf_file:
                    defaults_config.write(cf_file)
                return True
            elif answer in ["n", "no", "N", "No", "NO"]:
                print(
                    f"Update aborted, continuing with motif sequence version {defaults_config["VERSIONS"]['motifs']}"
                )
                if force_alg_update:
                    print("Update algorithms is set, updating algorithms...")
                    setup.update_algorithms()
                return True
            else:
                raise ValueError(
                    "Ineliglbe input, please answer the question with Yes/Y/yes/y or No/N/no/n."
                )
        else:
            print("RNA 3D motif sequences are up to date")
            if force_alg_update:
                print("Update algorithms is set, loading new motif sequences into algorithms...")
                setup.update_algorithms()
            return True
    except ValueError as error:
        print(error)
        return error


def _interactive_session(
    runtime_arguments: Namespace | configparser.ConfigParser,
) -> list[algorithm_output | None]:
    """Function is an infinite while Loop that always does one prediction, returns the results and waits for a new input."""
    result_list = []
    if isinstance(runtime_arguments, Namespace):
        proc_obj = proc.bgap_rna.from_argparse(runtime_arguments)
        output_file = runtime_arguments.output
        pool_boys = runtime_arguments.workers
        csv_separator = runtime_arguments.separator
    elif isinstance(runtime_arguments, configparser.ConfigParser):
        for argument in (x for x in runtime_arguments.default_section if x == ""):
            runtime_arguments.set(runtime_arguments.default_section, argument, None)
        proc_obj = proc.bgap_rna.from_config(runtime_arguments)
        output_file = runtime_arguments.get("VARIABLES", "output")
        pool_boys = runtime_arguments.getint("VARIABLES", "workers")
        csv_separator = runtime_arguments.get("VARIABLES", "separator")
    else:
        raise TypeError(
            "runtime args is neither a argparse Namespace nor a configparser.Configparser instance. Exiting..."
        )
    while True:
        print("Awaiting input...")
        user_input = input()
        if user_input in ["exit", "Exit", "Eixt", "Exi", "eixt"]:
            print("Exiting, thank you for using RNAmotiFold!")
            break
        else:
            proc_obj.input = user_input
            result = proc_obj.auto_run(
                o_file=output_file, pool_workers=pool_boys, output_csv_separator=csv_separator
            )
            result_list.append(result)
    return result_list  # Added result outputting just in case I wanna do something with that down the line.


def _uninteractive_session(
    runtime_arguments: Namespace | configparser.ConfigParser,
) -> list[algorithm_output]:
    if isinstance(runtime_arguments, Namespace):
        proc_obj = proc.bgap_rna.from_argparse(runtime_arguments)
        output_file = runtime_arguments.output
        pool_boys = runtime_arguments.workers
        csv_seaparator = runtime_arguments.separator
    elif isinstance(runtime_arguments, configparser.ConfigParser):
        for argument in (x for x in runtime_arguments.default_section if x == ""):
            runtime_arguments.set(runtime_arguments.default_section, argument, None)
        proc_obj = proc.bgap_rna.from_config(runtime_arguments)
        output_file = runtime_arguments.get("VARIABLES", "output")
        pool_boys = runtime_arguments.getint("VARIABLES", "workers")
        csv_seaparator = runtime_arguments.get("VARIABLES", "separator")
    result = proc_obj.auto_run(
        o_file=output_file,
        pool_workers=pool_boys,
        output_csv_separator=csv_seaparator,
    )
    return (
        result  # Added result outputting just in case I wanna do something with that down the line.
    )


def updates(no_update: bool, update_algorithms: bool):
    if not no_update:
        version_check_done = False
        while version_check_done is not True:
            try:
                version_check_done = _version_checks(update_algorithms)
            except ValueError as error:
                raise error
    else:
        print("Motif sequence updating disabled, continuing without update....")


def _duplicate_warning(mot_source: int, mot_orient: int) -> list[dict]:
    """Gives out a warning about possible duplicates in the currently used motif set"""
    if mot_source < 4:
        motif_set_in_use = get_current_set(mot_source, mot_orient)
        for motif_set in motif_set_in_use:
            for key in motif_set.keys():
                if key != "motif_mode":
                    print(
                        f"Warning possibles duplicates: {key}: {"".join(motif_set[key])}/{key.upper()}"
                    )
    else:
        print(
            "Duplicate warnings are disabled for custom motif sets, I trust you know what you're doing."
        )


# builds the names of the 3 motif sets in use currently from the -Q and -b input arguments
def get_current_set(mot_src: int, mot_ori: int) -> list[dict]:
    types = ["hairpins", "internals", "bulges"]
    match mot_src:
        case 1:
            src = "rna3d"
        case 2:
            src = "rfam"
        case 3:
            src = "both"
    match mot_ori:
        case 1:
            orient = "fw.csv"
        case 2:
            orient = "rv.csv"
        case 3:
            orient = "both.csv"
    filenames = ["_".join([src, x, orient]) for x in types]
    with open(
        Path(__file__).resolve().parent.joinpath("src", "data", "duplicates.json"), "r"
    ) as file:
        read_json = json.load(file, object_hook=dict)
    current_set = [x for x in read_json if x["motif_mode"] in filenames]
    return current_set


if __name__ == "__main__":
    print(
        "Starting RNAmotiFold, loading data from /src/data/config.ini and checking RNA 3D motif sequence versions ..."
    )
    rt_args = args.get_cmdarguments()
    try:
        if isinstance(rt_args, configparser.ConfigParser):
            updates(
                rt_args.getboolean("VARIABLES", "no_update"),
                rt_args.getboolean("VARIABLES", "update_algorithms"),
            )
            _duplicate_warning(
                rt_args.getint("VARIABLES", "motif_source"),
                rt_args.getint("VARIABLES", "motif_orientation"),
            )
            user_input = rt_args.get("VARIABLES", "input")
        elif isinstance(rt_args, Namespace):
            updates(rt_args.no_update, rt_args.update_algorithms)
            _duplicate_warning(rt_args.motif_source, rt_args.motif_orientation)
            user_input = rt_args.input
    except ValueError as error:
        raise error

    if len(user_input) > 0:
        out = _uninteractive_session(rt_args)

    else:
        out = _interactive_session(rt_args)
