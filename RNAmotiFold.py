from src.bgap_rna import bgap_rna, input_handler, motif_handler, subprocess_handler
from src.results import base_result,algorithm_output
import installer
import src.bgap_rna.alg_setup as setup
import src.input.arg_parsing as arg_parsing
import logging
from pathlib import Path
import sys
import copy

logger = logging.getLogger("RNAmotiFold")

try:
    import submodules.RNALoops.Misc.Applications.RNAmotiFold.motifs.get_RNA3D_motifs as motifs
except ImportError as e:
    raise ImportError(
        f"Submodule RNALoops was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {Path(__file__).absolute().parent}"
    )

def combine_calls(base_call: str, motif_subcalls: list[str]):
    if len(motif_subcalls) == 0:
        return [base_call]
    else:
        return [base_call + " " + x for x in motif_subcalls]


def check_install() -> bool:
    checkpath = Path(__file__).parent / "Build" / "bin" / "RNAmotiFold"
    return checkpath.is_file()

def create_inputs(calls:list[str],inputs:list[input_handler.algorithm_input]) -> list[input_handler.algorithm_input]:
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


if __name__ == "__main__":
    #Parse CMD and configure logger and result objects
    rt_args, additional_parameters = arg_parsing.get_cmdarguments()
    configure_logs(rt_args.loglevel,rt_args.logfile)
    base_result.result.separator = rt_args.separator
    logger.debug(rt_args)
    #Check if RNAmotiFold is installed and do updates if necessary/wanted
    if not check_install():
        installer.main()
        if not check_install():
            raise FileNotFoundError(
                "Could not find installed algorithm binaries, please run installer.py if you haven't yet"
            )
    else:
        if not rt_args.no_update:
            try:
                updated = setup.updates(motif_version=rt_args.version)
            except:
                pass
            else:
                if updated:
                    logger.debug(f"Updated to {rt_args.version}")
                else:
                    logger.debug(f"Failed to update to version {rt_args.version}")
        else:
            if rt_args.version == "current":
                rt_args.version = motifs.currently_installed()
    #Create all the support class instances to separately handle inputs, motif calls, algorithm calls and subprocesses
    input_maker = input_handler.input_handler(rt_args.process_type, rt_args.input)
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
        rt_args.workers,  rt_args.output, rt_args.alg_type(),len(motif_calls)
    )
    full_calls: list[str] = combine_calls(call_maker.call, motif_calls)

    if rt_args.input is not None:
        inputs: list[input_handler.algorithm_input] = input_maker.read_input(
            rt_args.process_type, rt_args.input, rt_args.id
        )
        alg_input = create_inputs(full_calls,inputs)
        results = subprocess_manager.run(alg_input,rt_args.fast_mode_merge)
    else:
        while True:
            print("Awaiting input...")
            user_input = input()
            if user_input.strip().lower() in ["exit","eixt","exi"]:
                print("Exiting...")
                break
            else:
                try:
                    rt_input = input_maker.read_input(rt_args.process_type,user_input,rt_args.id)
                except ValueError as e:
                    print(e)
                    continue
                except OSError as e:
                    print(e)
                    continue
                alg_input = create_inputs(full_calls,rt_input)
                results = subprocess_manager.run(alg_input,rt_args.fast_mode_merge)
    input_handler.algorithm_input.cleanup_temps()
    motif_subcall_maker.cleanup_tmp_files()