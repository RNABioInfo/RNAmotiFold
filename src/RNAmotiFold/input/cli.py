from importlib.util import module_from_spec, spec_from_file_location

from RNAmotiFold.bgap_rna import (
    bgap_rna,
    input_handler,
    motif_handler,
    subprocess_handler,
)
from RNAmotiFold.results import base_result, algorithm_output
import RNAmotiFold.bgap_rna.alg_setup as setup
import RNAmotiFold.input.arg_parsing as arg_parsing
import logging
from pathlib import Path
import sys
import copy
import subprocess
import os
import shutil

logger = logging.getLogger(__name__)

try:
    script_dir= setup.ROOT_DIR / "submodules" / "RNALoops" / "Misc" / "Applications" / "RNAmotiFold" / "motifs" / "get_RNA3D_motifs.py"
    spec = spec_from_file_location("uniteractive_update",script_dir)
    if spec is None or spec.loader is None:
        raise ImportError(f"Submodule RNALoops was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {setup.ROOT_DIR}")
    motifs = module_from_spec(spec)
    spec.loader.exec_module(motifs)
except ImportError as e:
    raise e

def combine_calls(base_call: str, motif_subcalls: list[str]):
    if len(motif_subcalls) == 0:
        return [base_call]
    else:
        return [base_call + " " + x for x in motif_subcalls]


def check_install() -> bool:
    checkpath = Path(__file__).parents[3].resolve() / "Build" / "bin" / "RNAmotiFold"
    return checkpath.exists()

def temp_cleanup():
    filepath = Path(__file__).resolve().parents[1]
    tempfolders = [x[0] for x in os.walk(filepath) if "tmp_" in x[0]]
    if len(tempfolders) > 0:
        logger.debug(f"Identified leftover tmp folder(s) from previous run: {", ".join(tempfolders)}, deleting...")
        for folder in tempfolders:
            try:
                shutil.rmtree(folder)
            except FileNotFoundError as e:
                logger.info(f"Could not delete some temp files from previous runs {folder}. Continuing without deleting it. This has no impact on the current run.")


def create_inputs(
    calls: list[str], inputs: list[input_handler.algorithm_input]
) -> list[input_handler.algorithm_input]:
    full_inputs: list[input_handler.algorithm_input] = []
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


def main() -> list[algorithm_output.algorithm_output | algorithm_output.error]:
    rt_args, additional_parameters = arg_parsing.get_cmdarguments()
    configure_logs(loglevel=rt_args.loglevel, logfile=rt_args.logfile)
    temp_cleanup()
    base_result.result.separator = rt_args.separator
    logger.debug(rt_args)
    # Check if RNAmotiFold is installed and do updates if necessary/wanted
    if not check_install():
        logger.critical(
            "Couldn't find RNAmotiFold, attempting to install algorithms and gapc if necessary"
        )
        setup.main(
            rt_args.version,
            rt_args.gapc_path,
            rt_args.perl_path,
            workers=rt_args.workers,
            cmake_path=rt_args.cmake_path,
        )
        if not check_install():
            raise FileNotFoundError(
                "Something went wrong setting up algorithms, check if dependencies are installed and re-run installer.py"
            )
        else:
            logger.critical("Installation successfull, running RNAmotiFold")
    else:
        if rt_args.update:
            logger.info(
                "Update is set, attempting to update algorithms to given version or current version"
            )
            try:
                updated: bool = setup.updates(motif_version=rt_args.version)
            except RuntimeError as e:
                raise e
            except subprocess.CalledProcessError as e:
                raise e
            else:
                if updated:
                    logger.info(f"Updated to {rt_args.version}")
                else:
                    logger.info(
                        f"Failed to update to version {rt_args.version}, trying to run with currently installed version"
                    )
        else:
            if rt_args.version == "current":
                rt_args.version = motifs.currently_installed()
            else:
                if rt_args.version != motifs.currently_installed():
                    updated: bool = setup.updates(motif_version=rt_args.version)
    # Create all the support class instances to separately handle inputs, motif calls, algorithm calls and subprocesses
    input_maker = input_handler.input_handler(
        process_type=rt_args.process_type, user_input=rt_args.input
    )
    motif_subcall_maker = motif_handler.motif_handler(
        rt_args.motif_list,
        rt_args.fast_mode,
        rt_args.custom_hairpins,
        rt_args.custom_internals,
        rt_args.custom_bulges,
        rt_args.replace_hairpins,
        rt_args.replace_internals,
        rt_args.replace_bulges,
        rt_args.version,
    )
    motif_calls = motif_subcall_maker.motif_calls
    call_maker = bgap_rna.bgap_rna.from_script_parameters(rt_args)
    subprocess_manager = subprocess_handler.subprocess_handler(
        rt_args.workers,
        rt_args.output,
        rt_args.alg_type(),
        motif_subcall_maker.call_number,
    )
    full_calls: list[str] = combine_calls(call_maker.call, motif_calls)
    if rt_args.input is not None:
        inputs: list[input_handler.algorithm_input] = input_maker.read_input(
            process_type=rt_args.process_type,
            user_input=rt_args.input,
            id=rt_args.id,
        )
        alg_input = create_inputs(calls=full_calls, inputs=inputs)
        results = subprocess_manager.run(
            inputs=alg_input, merge_mfe_outputs=rt_args.fast_mode_merge
        )
    else:
        results:list[algorithm_output.algorithm_output|algorithm_output.error] = []
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
                alg_input: list[input_handler.algorithm_input] = create_inputs(full_calls, rt_input)
                results.extend(subprocess_manager.run(alg_input, rt_args.fast_mode_merge))
    input_handler.algorithm_input.cleanup_temps()
    motif_subcall_maker.cleanup_tmp_files()
    return results


if __name__ == "__main__":
    main()
