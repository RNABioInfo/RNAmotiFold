import logging
import os
import sys
from pathlib import Path
from typing import NoReturn

from RNAmotiFold.bgap_rna import alg_setup
from RNAmotiFold.bgap_rna.call_handler import main_call_from_script_params
from RNAmotiFold.bgap_rna.input_handler import (
    InputException,
    ReadyAlgorithmInput,
    get_input,
)
from RNAmotiFold.bgap_rna.motif_handler import (
    MotifCall,
    create_motif_subcalls,
)
from RNAmotiFold.bgap_rna.subprocess_handler import (
    make_handler,
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
        x
        for x in [entry.value for entry in AvailableBinaries]
        if not check_install(x)
    ]


def check_install(algorithm: str) -> bool:
    checkpath: Path = RNAMOTIFOLD_ROOT_DIR / "bin" / algorithm
    return checkpath.is_file() and os.access(checkpath, os.X_OK)


# configures all loggers with logging.basicConfig to use the same loglevel and output to the same destination
def _configure_logs(loglevel: str | int, logfile: Path | None) -> None:
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

    # Check if everything is set up properly
    not_installed = _check_all_algorithms()
    if len(not_installed) > 0:
        raise FileNotFoundError(
            f"Algorithms {','.join(not_installed)} is currently not installed, please run RNAmotiFold_installer"
        )

    # Create the base RNAmotiFold call and motif sub-calls
    main_call = main_call_from_script_params(rt_args)
    motif_calls: list[MotifCall] = create_motif_subcalls(
        rt_args.custom_hairpins,
        rt_args.custom_internals,
        rt_args.custom_bulges,
        rt_args.replace_hairpins,
        rt_args.replace_internals,
        rt_args.replace_bulges,
        rt_args.version,
        rt_args.motif_list,
        rt_args.single_motif_mode,
    )

    # Set up a subprocess Manager
    subprocess_manager = make_handler(
        rt_args.workers, rt_args.output, rt_args.alg_type
    )

    if rt_args.input is not None:
        ready_inputs: list[ReadyAlgorithmInput] = get_input(
            main_call,
            motif_calls,
            rt_args.process_type,
            rt_args.input,
            rt_args.id,
        )
        result = subprocess_manager.auto_run(
            inputs=ready_inputs,
            no_print=False,
        )
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
                    rt_input: list[ReadyAlgorithmInput] = get_input(
                        main_call,
                        motif_calls,
                        rt_args.process_type,
                        user_input,
                        rt_args.id,
                    )
                except ValueError as e:
                    print(e)
                    continue
                except OSError as e:
                    print(e)
                    continue
                except InputException as e:
                    print(e)
                    continue
                results.append(
                    subprocess_manager.auto_run(
                        inputs=rt_input,
                        no_print=False,
                    )
                )
        result = AlgorithmOutput.merge_list(results)
    TempFileHandler.cleanup()
    return 1 if result.errors else 0


def _rnamotifold_install() -> NoReturn:
    args = arg_parsing.get_install_arguments()
    _configure_logs("DEBUG", None)
    not_installed = _check_all_algorithms()
    if args.force or len(not_installed) > 0:
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
    _rnamotifold_install()
