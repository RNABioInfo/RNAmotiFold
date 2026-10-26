import copy
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Literal, NoReturn

from RNAmotiFold.bgap_rna import alg_setup
from RNAmotiFold.bgap_rna.call_handler import CallHandler
from RNAmotiFold.bgap_rna.input_handler import (
    AlgorithmInput,
    InputHandler,
)
from RNAmotiFold.bgap_rna.motif_handler import MotifHandler
from RNAmotiFold.bgap_rna.subprocess_handler import SubprocessHandler
from RNAmotiFold.input import action_overwrites, arg_parsing
from RNAmotiFold.results.outputs import (
    AlgorithmOutput,
    ErrorOutput,
)
from RNAmotiFold.results.results import BaseResult

from RNAmotiFold import AVAILABLE_BINARIES, RNAMOTIFOLD_ROOT_DIR

logger = logging.getLogger(__name__)


def _check_all_algorithms() -> list[str]:
    """Checks if all algorithms are installed and returns True if they are, False if not"""
    return [x for x in AVAILABLE_BINARIES if not check_install(x)]


def combine_calls(base_call: str, motif_subcalls: list[str]) -> list[str]:
    if len(motif_subcalls) == 0:
        return [base_call]
    else:
        return [base_call + " " + x for x in motif_subcalls]


def check_install(algorithm: str) -> bool:
    checkpath = RNAMOTIFOLD_ROOT_DIR / "bin" / algorithm
    return checkpath.exists()


def _temp_cleanup() -> None:
    filepath = Path(__file__).resolve().parents[1]
    tempfolders = [x[0] for x in os.walk(filepath) if "tmp_" in x[0]]
    if len(tempfolders) > 0:
        logger.debug(
            f"Identified leftover tmp folder(s) from previous run: {', '.join(tempfolders)}, deleting..."
        )
        for folder in tempfolders:
            try:
                shutil.rmtree(folder)
            except FileNotFoundError:
                logger.info(
                    f"Could not delete some temp files from previous runs {folder}. Continuing without deleting it. This has no impact on the current run."
                )


def create_inputs(
    calls: list[str], inputs: list[AlgorithmInput]
) -> list[AlgorithmInput]:
    full_inputs: list[AlgorithmInput] = []
    for input in inputs:
        for call in calls:
            new_input = copy.copy(input)
            new_input.call = call
            full_inputs.append(new_input)
    return full_inputs


# configures all loggers with logging.basicConfig to use the same loglevel and output to the same destination
def _configure_logs(loglevel: str, logfile: Path | None) -> None:
    if logfile is not None:
        logging.basicConfig(
            filename=logfile,
            filemode="a+",
            level=loglevel,
            format="%(asctime)s:%(name)s:%(levelname)s:%(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    else:
        logging.basicConfig(
            stream=sys.stderr,
            level=loglevel,
            format="%(asctime)s:%(name)s:%(levelname)s:%(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )


def _main(
    algorithm: (
        Literal["rnamotifold", "rnamotices", "rnamotialign", "rnamosh"]
        | None
    ) = None,
) -> int:
    rt_args: arg_parsing.ScriptParameters = arg_parsing.get_cmdarguments()
    if algorithm is not None:
        rt_args.algorithm = (
            action_overwrites.AlgorithmMatching.algorithm_matching_function(
                algorithm
            )
        )  # Corrects capitalization
    _configure_logs(loglevel=rt_args.loglevel, logfile=rt_args.logfile)
    _temp_cleanup()
    BaseResult.separator = rt_args.separator
    logger.debug(rt_args)

    # Create all the support class instances to separately handle inputs, motif calls, algorithm calls and subprocesses
    input_maker = InputHandler(
        process_type=rt_args.process_type, user_input=rt_args.input
    )
    motif_subcall_maker = MotifHandler.from_script_parameters(rt_args)
    call_maker = CallHandler.from_script_parameters(rt_args)
    subprocess_manager = SubprocessHandler.from_script_parameters(rt_args)

    # Combine the (potentially) separate motif calls and the general commandline call
    full_calls: list[str] = combine_calls(
        call_maker.call, motif_subcall_maker.motif_calls
    )
    not_installed = _check_all_algorithms()

    # Check if RNAmotiFold is installed and do updates if necessary/wanted
    if call_maker.algorithm_binary in not_installed:
        raise FileNotFoundError(
            f"Algorithm {call_maker.algorithm_binary} is currently not installed, please re-run RNAmotiFold_installer"
        )

    if rt_args.input is not None:
        inputs: list[AlgorithmInput] = input_maker.read_input(
            process_type=rt_args.process_type,
            user_input=rt_args.input,
            id=rt_args.id,
        )
        cmd_inputs: list[AlgorithmInput] = create_inputs(
            calls=full_calls, inputs=inputs
        )
        results = [
            subprocess_manager.auto_run(
                inputs=cmd_inputs,
                calls_per_input=motif_subcall_maker.call_number,
                no_print=False,
            )
        ]
    else:
        results: list[AlgorithmOutput] = []
        while True:
            print("Awaiting input...")
            user_input = input()
            if user_input.strip().lower() in ["exit", "eixt", "exi"]:
                print("Exiting...")
                break
            else:
                try:
                    rt_input = input_maker.read_input(
                        rt_args.process_type, user_input, rt_args.id
                    )
                except ValueError as e:
                    print(e)
                    continue
                except OSError as e:
                    print(e)
                    continue
                alg_input: list[AlgorithmInput] = create_inputs(
                    full_calls, rt_input
                )
                results.append(
                    subprocess_manager.auto_run(
                        inputs=alg_input,
                        calls_per_input=motif_subcall_maker.call_number,
                        no_print=False,
                    )
                )
    AlgorithmInput.cleanup_temps()
    motif_subcall_maker.cleanup_tmp_files()
    if any(isinstance(x, ErrorOutput) for x in results):
        return 1
    return 0


def _rnamotifold_install() -> NoReturn:
    args = arg_parsing.get_install_arguments()
    not_installed = _check_all_algorithms()
    if len(not_installed) == 0:
        if args.force:
            exit_code = alg_setup.main(
                AVAILABLE_BINARIES,
                args.gapcM_path,
                args.perl_path,
                args.make_path,
                args.workers,
            )
        else:
            print("Algorithms are already installed")
            exit_code = 0
    else:
        exit_code = alg_setup.main(
            AVAILABLE_BINARIES,
            args.gapcM_path,
            args.perl_path,
            args.make_path,
            args.workers,
        )
    sys.exit(exit_code)


def _rnamotifold_cli() -> NoReturn:
    """Main function to run RNAmotiFold from command line"""
    exit_code = _main("rnamotifold")
    sys.exit(exit_code)


def _rnamotices_cli() -> NoReturn:
    """Main function to run RNAmotiCes from command line"""
    exit_code = _main("rnamotices")
    sys.exit(exit_code)


def _rnamotialign_cli() -> NoReturn:
    """Main function to run RNAmotiAlign from command line"""
    exit_code = _main("rnamotialign")
    sys.exit(exit_code)


def _rnamosh_cli() -> NoReturn:
    """Main function to run RNAmoSh from command line"""
    exit_code = _main("rnamosh")
    sys.exit(exit_code)


if __name__ == "__main__":
    _main()
