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
from time import sleep
from multiprocessing import Process


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


def update_wheel():
    sleep(0.1)  # Prevents race condition on printing "Updating..."
    for frame in r"-\|/-\|/":
        print("\r", frame, sep="", end="", flush=True)
        sleep(0.2)


def _version_checks(force_alg_update: bool):
    defaults_config = configparser.ConfigParser(allow_no_value=True)
    defaults_config.read_file(open(args.script_parameters.defaults_config_path))
    try:
        update_needed = proc.bgap_rna.check_motif_versions(defaults_config["VERSIONS"]["motifs"])
        if update_needed:
            print(
                "There is a new set of RNA 3D Motif sequences available. You may need to update the motif.json file manually. Update RNAmotiFold ? [y/n]",
                end=" ",
            )
            answer = input()
            if answer.lower() in ["y", "ye", "yes"]:
                p = Process(target=setup.update_sequences_algorithms)
                p.start()
                while True:
                    update_wheel()
                    if not p.is_alive():
                        break
                p.join()
                defaults_config.set("VERSIONS", "motifs", proc.bgap_rna.get_current_motif_version())
                with open(args.script_parameters.defaults_config_path, "w") as file:
                    defaults_config.write(file)
                return True

            elif answer.lower() in ["n", "no"]:
                print(
                    f"Skipping update, continuing with motif sequence version {defaults_config['VERSIONS']['motifs']}"
                )
                if force_alg_update:
                    print("Update algorithms is set, updating algorithms...")  # make log
                    setup.update_algorithms()
                return True
            else:
                raise ValueError("Please answer the question with yes or no")
        else:
            print("RNA 3D motif sequences are up to date")  # make log
            if force_alg_update:
                print(
                    "Update algorithms is set, loading new motif sequences into algorithms..."
                )  # make log
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
        proc_obj = proc.bgap_rna.from_config(runtime_arguments, "VARIABLES")
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
        if user_input.strip() in ["exit", "Exit", "Eixt", "Exi", "eixt"]:
            print("Exiting, thank you for using RNAmotiFold!")
            break
        elif user_input in ["status", "Status", "state"]:
            print(proc_obj)
        elif user_input in ["h", "help", "-h"]:
            print(
                f"You are currently using the following algorithm call:\n{str(proc_obj)}\n Please input a RNA/DNA sequence or a fasta,fastq or stockholm formatted sequence file."
            )
        else:
            proc_obj.input = user_input.strip()
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
        output_file = runtime_arguments.output  # type:str
        pool_boys = runtime_arguments.workers  # type: int
        csv_seaparator = runtime_arguments.separator  # type:str
    elif isinstance(runtime_arguments, configparser.ConfigParser):
        proc_obj = proc.bgap_rna.from_config(runtime_arguments)
        output_file = runtime_arguments.get("VARIABLES", "output")  # type:str
        pool_boys = runtime_arguments.getint("VARIABLES", "workers")
        csv_seaparator = runtime_arguments.get("VARIABLES", "separator")  # type:str
    result = proc_obj.auto_run(
        o_file=output_file,
        pool_workers=pool_boys,
        output_csv_separator=csv_seaparator,
    )
    # Added result outputting just in case I wanna do something with that down the line.
    return result


def updates(no_update: bool, update_algorithms: bool):
    if not no_update:
        version_check_done = False
        while version_check_done is not True:
            try:
                version_check_done = _version_checks(update_algorithms)
            except ValueError as error:
                raise error
    else:
        print(
            "Motif sequence updating disabled, continuing without update..."
        )  # make into a log entry


def _duplicate_warning(mot_source: int, mot_orient: int) -> list[dict]:
    """Gives out a warning about possible duplicates in the currently used motif set"""
    if mot_source < 4:
        motif_set_in_use = get_current_set(mot_source, mot_orient)
        for motif_set in motif_set_in_use:
            for key in motif_set.keys():
                if len(key) == 1:
                    print(
                        f"Warning possibles duplicates: {key}: {''.join(motif_set[key])}/{key.upper()}"
                    )  # Also log this for uninteractive sessions.
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

    if user_input is not None:
        out = _uninteractive_session(rt_args)

    else:
        out = _interactive_session(rt_args)
