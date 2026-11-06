import logging
import os
import sys
from pathlib import Path
from typing import NoReturn

from RNAmotiFold.bgap_rna import alg_setup
from RNAmotiFold.bgap_rna.call_handler import CallHandler
from RNAmotiFold.bgap_rna.input_handler import (
    AlgorithmInput,
    ReadyAlgorithmInput,
    create_input_handler,
)
from RNAmotiFold.bgap_rna.motif_handler import MotifCall, MotifHandler
from RNAmotiFold.bgap_rna.subprocess_handler import (
    create_subprocess_handler,
)
from RNAmotiFold.bgap_rna.tempfile_handler import TempFileHandler
from RNAmotiFold.input import arg_parsing
from RNAmotiFold.Misc.EnumTypes import AlgorithmType, AvailableBinaries
from RNAmotiFold.results.outputs import AlgorithmOutput
from RNAmotiFold.results.results import BaseResult

from RNAmotiFold import RNAMOTIFOLD_ROOT_DIR

logger = logging.getLogger(__name__)


def _check_all_algorithms() -> list[str]:
    """Returns algoriths that are not installed"""
    return [
        x for x in [entry.value for entry in AvailableBinaries] if not check_install(x)
    ]


def check_install(algorithm: str) -> bool:
    checkpath: Path = RNAMOTIFOLD_ROOT_DIR / "bin" / algorithm
    return checkpath.is_file() and os.access(checkpath, os.X_OK)


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
    override_algorithm: AlgorithmType | None = None,
) -> int:
    rt_args: arg_parsing.ScriptParameters = arg_parsing.get_cmdarguments()
    if override_algorithm is not None:
        rt_args.algorithm = override_algorithm

    _configure_logs(loglevel=rt_args.loglevel, logfile=rt_args.logfile)

    BaseResult.separator = rt_args.separator
    logger.debug(rt_args)

    # Create all the support class instances to separately handle inputs, motif calls, algorithm calls and subprocesses
    input_maker = create_input_handler(process_type=rt_args.process_type)
    motif_subcall_maker = MotifHandler.from_script_parameters(rt_args)
    call_maker = CallHandler.from_script_parameters(rt_args)
    subprocess_manager = create_subprocess_handler(rt_args)

    # Create the necessary calls
    main_call = call_maker.make_main_call()
    motif_calls: list[MotifCall] = motif_subcall_maker.create_motif_calls()

    # Check if RNAmotiFold is installed and do updates if necessary/wanted
    not_installed = _check_all_algorithms()
    if call_maker.algorithm_binary in not_installed:
        raise FileNotFoundError(
            f"Algorithm {call_maker.algorithm_binary} is currently not installed, please re-run RNAmotiFold_installer"
        )

    if rt_args.input is not None:
        inputs: list[AlgorithmInput] = input_maker.read_input(rt_args.input, rt_args.id)
        ready_inputs = [
            ReadyAlgorithmInput.create_ready_input(x, main_call, motif_calls)
            for x in inputs
        ]
        result = subprocess_manager.auto_run(
            inputs=ready_inputs,
            no_print=False,
        )
        return 1 if result.errors else 0
    else:
        results: list[AlgorithmOutput] = []
        while True:
            try:
                print("Awaiting input...")
                user_input = input()
            except EOFError:
                print(
                    "Exiting due to EOF Error when attempting to get input from the commandline"
                )
                break
            if user_input.strip().casefold() == "exit":
                break
            else:
                try:
                    rt_input: list[AlgorithmInput] = input_maker.read_input(
                        user_input, rt_args.id
                    )
                except ValueError as e:
                    print(e)
                    continue
                except OSError as e:
                    print(e)
                    continue
                alg_input: list[ReadyAlgorithmInput] = [
                    ReadyAlgorithmInput.create_ready_input(x, main_call, motif_calls)
                    for x in rt_input
                ]
                results.append(
                    subprocess_manager.auto_run(
                        inputs=alg_input,
                        no_print=False,
                    )
                )
    result = AlgorithmOutput.merge_list(results)
    TempFileHandler.cleanup()
    return 1 if result.errors else 0


def _rnamotifold_install() -> NoReturn:
    args = arg_parsing.get_install_arguments()
    not_installed = _check_all_algorithms()
    if len(not_installed) == 0:
        if args.force:
            exit_code = alg_setup.main(
                [entry.value for entry in AvailableBinaries],
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
            [entry.value for entry in AvailableBinaries],
            args.gapcM_path,
            args.perl_path,
            args.make_path,
            args.workers,
        )
    sys.exit(exit_code)


def _rnamotifold_cli() -> NoReturn:
    """Main function to run RNAmotiFold from command line"""
    exit_code = _main(AlgorithmType.RNAMOTIFOLD)
    sys.exit(exit_code)


def _rnamotices_cli() -> NoReturn:
    """Main function to run RNAmotiCes from command line"""
    exit_code = _main(AlgorithmType.RNAMOTICES)
    sys.exit(exit_code)


def _rnamotialign_cli() -> NoReturn:
    """Main function to run RNAmotiAlign from command line"""
    exit_code = _main(AlgorithmType.RNAMOTIALIGN)
    sys.exit(exit_code)


def _rnamosh_cli() -> NoReturn:
    """Main function to run RNAmoSh from command line"""
    exit_code = _main(AlgorithmType.RNAMOSH)
    sys.exit(exit_code)


if __name__ == "__main__":
    _main()
