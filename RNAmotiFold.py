from src.bgap_rna import bgap_rna, input_handler, motif_handler, subprocess_handler
import installer
import src.bgap_rna.alg_setup as setup
import src.input.arg_parsing as arg_parsing
import logging
from pathlib import Path
import sys

logger = logging.getLogger("RNAmotiFold")

try:
    import submodules.RNALoops.Misc.Applications.RNAmotiFold.motifs.get_RNA3D_motifs as motifs
except ImportError as e:
    raise ImportError(
        f"Submodule RNALoops was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {Path(__file__).absolute().parent}"
    )


# Interactive session to run multiple predictions in an "interactive" environment
# def _interactive_session(
#    runtime_arguments: ScriptParameters,
# ) -> list[algorithm_output | error]:
#    """Function is an infinite while Loop that always does one prediction, appends the result to a list and waits for a new input. List of results is returned"""
#    result_list: list[list[algorithm_output | error]] = []
#    proc_obj = bgap_rna.bgap_rna.from_script_parameters(runtime_arguments)
#    logger.debug("Created bgap_rna obj: " + repr(proc_obj))
#    while True:
#        print("Awaiting input...")
#        user_input = input()
#        if user_input.strip().lower() in ["exit", "eixt", "exi"]:
#            logger.debug("Exit was given as input, exiting...")
#            break
#        elif user_input.strip().lower() in ["h", "help", "-h"]:
#            print(
#                f"You are currently using the following algorithm call:\n{str(proc_obj)}\n Please input a RNA/DNA sequence or a fasta, fastq or stockholm formatted sequence file."
#            )
#        else:
#            try:
#                realtime_input: (
#                    FastaIO.FastaIterator
#                    | QualityIO.FastqPhredIterator
#                    | Generator[SeqRecord, None, None]
#                    | SeqRecord
#                    | list[SeqRecord]
#                ) = _input_check(user_input, runtime_arguments.id)
#            except ValueError as v_error:
#                print(v_error)
#            except OSError as os_error:
#                print(os_error)
#            #else:
#    result: list[
#        algorithm_output | error
#    ] = proc_obj.auto_run(
#        realtime_input,
#        version=runtime_arguments.version,
#        o_file=runtime_arguments.output,
#        pool_workers=runtime_arguments.workers,
#        output_csv_separator=runtime_arguments.separator,
#        merge=runtime_arguments.fast_mode_merge,
#    )
#    result_list.append(result)
# flat_list = results.flatten(result_list)
# proc_obj.cleanup_temp_files()
# return flat_list  # Added result outputting just in case I wanna do something with that down the line.


# Uninteractive session in case of preset input, just does the calculation and exits
# def _uninteractive_session(
#    runtime_arguments: arg_parsing.script_parameters,
# ) -> list[results.algorithm_output | results.error]:
#    runtime_input = _input_check(runtime_arguments.input, runtime_arguments.id)  # type: ignore cause we can only get here by argument not being None in main
#    proc_obj = bgap.bgap_rna.from_script_parameters(runtime_arguments)
#    logger.debug("Created bgap_rna obj: " + repr(proc_obj))
#    result: list[results.algorithm_output | results.error] = (
#        proc_obj.auto_run(
#            user_input=runtime_input,
#            version=runtime_arguments.version,
#            o_file=runtime_arguments.output,
#            pool_workers=runtime_arguments.workers,
#            output_csv_separator=runtime_arguments.separator,
#            merge=runtime_arguments.fast_mode_merge,
#            name=runtime_arguments.id,
#        )
#    )
#    proc_obj.cleanup_temp_files()
#    return result


def combine_calls(base_call: str, motif_subcalls: list[str]):
    if len(motif_subcalls) == 0:
        return [base_call]
    else:
        return [base_call + " " + x for x in motif_subcalls]


def check_install() -> bool:
    checkpath = Path(__file__).parent / "Build" / "bin" / "RNAmotiFold"
    return checkpath.is_file()


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
    rt_args, additional_parameters = arg_parsing.get_cmdarguments()
    if not check_install():
        installer.main()
        if not check_install():
            raise FileNotFoundError(
                "Could not find installed algorithm binaries, please run installer.py if you haven't yet"
            )
    else:
        if not rt_args.no_update:
            setup.updates(motif_version=rt_args.version)

    rt_args.version = motifs.currently_installed().replace(".", "_")
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
    call_maker = bgap_rna.bgap_rna.from_script_parameters(rt_args)
    subprocess_manager = subprocess_handler.subprocess_handler(
        rt_args.workers, None, rt_args.output, rt_args.alg_type()
    )
    if rt_args.input is not None:
        inputs: list[input_handler.algorithm_input] = input_maker.read_input(
            rt_args.process_type, rt_args.input, rt_args.id
        )
    else:
        raise ValueError("no input set")
    full_calls: list[str] = combine_calls(call_maker.call, motif_subcall_maker.motif_calls)
    full_inputs: list[input_handler.algorithm_input] = []
    for call in full_calls:
        for input in inputs:
            input.call = call
            full_inputs.append(input)
    subprocess_manager.inputs = full_inputs
    results = subprocess_manager.run()

    input_handler.algorithm_input.cleanup_temps()
    motif_subcall_maker.cleanup_tmp_files()
    # try:
    #    configure_logs(
    #        loglevel=rt_args.loglevel, logfile=rt_args.logfile
    #    )
    #    if not rt_args.no_update:
    #        setup.updates(motif_version=rt_args.version)
    #    rt_args.version = motifs.currently_installed().replace(".", "_")
    # except ValueError as error:
    #    raise error
    # logger.debug("Input args: " + repr(rt_args))
    # if rt_args.input is not None:
    #    logger.info("Input is set, starting calculations")
    #    out: list[results.algorithm_output | results.error] = (
    #        _uninteractive_session(runtime_arguments=rt_args)
    #    )
    # else:
    #    logger.info("No input set, starting interactive session")
    #    out: list[results.algorithm_output | results.error] = (
    #        _interactive_session(runtime_arguments=rt_args)
    #    )
