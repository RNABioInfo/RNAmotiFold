from RNAmotiFold.input import action_overwrites
from RNAmotiFold.results.base_result import Result
from RNAmotiFold.results.algorithm_output import AlgorithmOutput, AlgorithmError
import RNAmotiFold.bgap_rna.alg_setup as alg_setup
import RNAmotiFold.input.arg_parsing as arg_parsing
from RNAmotiFold.bgap_rna.input_handler import InputHandler, AlgorithmInput
from RNAmotiFold.bgap_rna.motif_handler import MotifHandler
from RNAmotiFold.bgap_rna.call_handler import CallHandler
from RNAmotiFold.bgap_rna.subprocess_handler import SubprocessHandler
from RNAmotiFold import AVAILABLE_ALGORITHMS
import RNAmotiFold
import logging
from pathlib import Path
import sys
from typing import Literal
import copy
import os
import shutil

logger = logging.getLogger(__name__)


def check_all_algorithms() -> list[str]:
    """Checks if all algorithms are installed and returns True if they are, False if not"""
    return [x for x in AVAILABLE_ALGORITHMS if not check_install(x)]


def combine_calls(base_call: str, motif_subcalls: list[str]):
    if len(motif_subcalls) == 0:
        return [base_call]
    else:
        return [base_call + " " + x for x in motif_subcalls]


def check_install(algorithm: str) -> bool:
    checkpath = RNAmotiFold.RNAMOTIFOLD_ROOT_DIR / "bin" / algorithm
    return checkpath.exists()


def temp_cleanup():
    filepath = Path(__file__).resolve().parents[1]
    tempfolders = [x[0] for x in os.walk(filepath) if "tmp_" in x[0]]
    if len(tempfolders) > 0:
        logger.debug(
            f"Identified leftover tmp folder(s) from previous run: {", ".join(tempfolders)}, deleting..."
        )
        for folder in tempfolders:
            try:
                shutil.rmtree(folder)
            except FileNotFoundError as e:
                logger.info(
                    f"Could not delete some temp files from previous runs {folder}. Continuing without deleting it. This has no impact on the current run."
                )


def create_inputs(calls: list[str], inputs: list[AlgorithmInput]) -> list[AlgorithmInput]:
    full_inputs: list[AlgorithmInput] = []
    for input in inputs:
        for call in calls:
            new_input = copy.copy(input)
            new_input.call = call
            full_inputs.append(new_input)
    return full_inputs


# configures all loggers with logging.basicConfig to use the same loglevel and output to the same destination
def configure_logs(loglevel: str, logfile: Path | None) -> None:
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


def main(
    algorithm: Literal["rnamotifold", "rnamotices", "rnamotialign", "rnamosh"] | None = None,
) -> int:
    rt_args: arg_parsing.ScriptParameters = arg_parsing.get_cmdarguments()
    if algorithm is not None:
        rt_args.algorithm = action_overwrites.AlgorithmMatching.algorithm_matching_function(
            algorithm
        )  # Corrects capitalization
    configure_logs(loglevel=rt_args.loglevel, logfile=rt_args.logfile)
    temp_cleanup()
    Result.separator = rt_args.separator
    logger.debug(rt_args)
    # Create all the support class instances to separately handle inputs, motif calls, algorithm calls and subprocesses

    input_maker = InputHandler(process_type=rt_args.process_type, user_input=rt_args.input)
    motif_subcall_maker = MotifHandler.from_script_parameters(rt_args)

    call_maker = CallHandler.from_script_parameters(rt_args)
    subprocess_manager = SubprocessHandler.from_script_parameters(rt_args)
    full_calls: list[str] = combine_calls(call_maker.call, motif_subcall_maker.motif_calls)

    not_installed = check_all_algorithms()

    # Check if RNAmotiFold is installed and do updates if necessary/wanted
    if len(not_installed) > 0:
        logger.critical(f"Compiling algorithms {[", ".join(not_installed)]}...")
        alg_setup.main(
            not_installed,
            rt_args.gapc_path,
            rt_args.perl_path,
            rt_args.workers,
        )
        if not check_install(call_maker.algorithm_binary):
            raise FileNotFoundError(
                f"Something went wrong setting up {call_maker.algorithm_binary}, check if dependencies are installed and re-run installer.py"
            )

    if rt_args.input is not None:
        inputs: list[AlgorithmInput] = input_maker.read_input(
            process_type=rt_args.process_type,
            user_input=rt_args.input,
            id=rt_args.id,
        )
        cmd_inputs: list[AlgorithmInput] = create_inputs(calls=full_calls, inputs=inputs)
        results = subprocess_manager.run(
            inputs=cmd_inputs,
            merge_mfe_outputs=rt_args.fast_mode_merge,
            calls_per_input=motif_subcall_maker.call_number,
            no_print=False,
        )
    else:
        results: list[AlgorithmOutput | AlgorithmError] = []
        while True:
            print("Awaiting input...")
            user_input = input()
            if user_input.strip().lower() in ["exit", "eixt", "exi"]:
                print("Exiting...")
                break
            else:
                try:
                    rt_input = input_maker.read_input(rt_args.process_type, user_input, rt_args.id)
                except ValueError as e:
                    print(e)
                    continue
                except OSError as e:
                    print(e)
                    continue
                alg_input: list[AlgorithmInput] = create_inputs(full_calls, rt_input)
                results.extend(subprocess_manager.run(alg_input, rt_args.fast_mode_merge))
    AlgorithmInput.cleanup_temps()
    motif_subcall_maker.cleanup_tmp_files()
    if all([isinstance(x, AlgorithmOutput) for x in results]):
        return 0
    return 1


def rnamotifold():
    """Main function to run RNAmotiFold from command line"""
    exit_code = main("rnamotifold")
    sys.exit(exit_code)


def rnamotices():
    """Main function to run RNAmotiCes from command line"""
    exit_code = main("rnamotices")
    sys.exit(exit_code)


def rnamotialign():
    """Main function to run RNAmotiAlign from command line"""
    exit_code = main("rnamotialign")
    sys.exit(exit_code)


def rnamosh():
    """Main function to run RNAmoSh from command line"""
    exit_code = main("rnamosh")
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
