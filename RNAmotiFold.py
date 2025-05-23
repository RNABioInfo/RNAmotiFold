import src.bgap_rna as bgap
import src.args as args
import src.results as results
import setup
import logging
from pathlib import Path
from typing import Generator
from Bio import SeqIO
from Bio.SeqIO import FastaIO, QualityIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import sys
from typing import Optional

logger = logging.getLogger("RNAmotiFold")


# Interactive session to run multiple predictions in an "interactive" environment
def _interactive_session(
    runtime_arguments: args.script_parameters,
) -> list[results.algorithm_output | results.error | list[results.algorithm_output]]:
    """Function is an infinite while Loop that always does one prediction, returns the results and waits for a new input."""
    result_list:list[results.algorithm_output | results.error | list[results.algorithm_output]] = []
    proc_obj = bgap.bgap_rna.from_script_parameters(runtime_arguments)
    logger.debug("Created bgap_rna obj: " + repr(proc_obj))
    while True:
        print("Awaiting input...")
        user_input = input()
        if user_input.strip().lower() in ["exit", "eixt", "exi"]:
            logger.debug("Exit was given as input, exiting...")
            break
        elif user_input.strip().lower() in ["h", "help", "-h"]:
            print(
                f"You are currently using the following algorithm call:\n{str(proc_obj)}\n Please input a RNA/DNA sequence or a fasta, fastq or stockholm formatted sequence file."
            )
        else:
            try:
                realtime_input = _input_check(user_input, runtime_arguments.id)
            except ValueError as error:
                print(error)
            else:
                result = proc_obj.auto_run(
                    realtime_input,
                    o_file=runtime_arguments.output,
                    pool_workers=runtime_arguments.workers,
                    output_csv_separator=runtime_arguments.separator,
                )
                result_list.append(result)
    return result_list  # Added result outputting just in case I wanna do something with that down the line.


# Uninteractive session in case of preset input, just does the calculation and exits
def _uninteractive_session(
    runtime_arguments: args.script_parameters,
) -> results.algorithm_output | results.error | list[results.algorithm_output]:
    runtime_input = _input_check(runtime_arguments.input, runtime_arguments.id) #type:ignore cause we can only get here by argument not being None in main
    proc_obj = bgap.bgap_rna.from_script_parameters(runtime_arguments)
    logger.debug("Created bgap_rna obj: " + repr(proc_obj))
    result = proc_obj.auto_run(
        user_input=runtime_input,
        o_file=runtime_arguments.output,
        pool_workers=runtime_arguments.workers,
        output_csv_separator=runtime_arguments.separator,
    )
    return result


# Finds File type based on file ending
def _find_filetype(file_path: Path) -> tuple[bool, str]:
    if file_path.suffixes[-1] == ".gz" or file_path.suffixes[-1] == ".zip":
        file_extension = file_path.suffixes[-2]
        input_zipped = True
    else:
        file_extension = file_path.suffixes[-1]
        input_zipped = False
    match file_extension:
        case (
            ".fasta"
            | ".fas"
            | ".fa"
            | ".fna"
            | ".ffn"
            | ".faa"
            | ".mpfa"
            | ".frn"
            | ".txt"
            | ".fsa"
        ):
            filetype = "fasta"

        case ".fastq" | ".fq":
            filetype = "fastq"
        case ".stk" | ".stockholm" | ".sto":
            filetype = "stockholm"
        case _:
            logger.critical(
                "Filetype was not recognized as fasta, fastq or stockholm format. Or file could not be unpacked, please ensure it is zipped with either .gz or .zip or not zipped at all."
            )
            raise TypeError("Unable to read given input file.")
    logger.debug(f"Recognized filetype as {filetype}.")
    return (input_zipped, filetype)


# Read input file
def _read_input_file(
    file_path: Path,
) -> (
    FastaIO.FastaIterator
    | QualityIO.FastqPhredIterator
    | Generator[SeqRecord, None, None]
):
    (zipped, filetype) = _find_filetype(file_path)
    if not zipped:
        return SeqIO.parse(file_path, filetype) #type:ignore
    else:
        with gzip.open(file_path, "rt") as handle:
            return SeqIO.parse(handle, filetype) #type:ignore


# This function still has a lot of leftover functionality from when it was part of the bgap_rna class, shouldn't really matter and I'll leave it in case I need it again later I guess.
def _input_check(user_input: str, id: str):
    if Path(user_input).resolve().is_file():
        logger.info("Recognized input as filepath, reading...")
        return _read_input_file(Path(user_input).resolve())
    else:
        if any(c not in "AUCGTaucgt+" for c in set(user_input)):
            raise ValueError(
                "Input string was neither a viable file path nor a viable RNA or DNA sequence"
            )
        else:
            return SeqRecord(seq=Seq(user_input), id=id)


# configures all loggers with logging.basicConfig to use the same loglevel and output to the same destination
def configure_logs(loglevel: str, logfile: Optional[Path]):
    if logfile is not None:
        logging.basicConfig(
            filename=logfile,
            filemode="a+",
            level=loglevel,
            format="%(asctime)s:%(name)s:%(levelname)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    else:
        logging.basicConfig(
            stream=sys.stderr,
            level=loglevel,
            format="%(asctime)s:%(name)s:%(levelname)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )


if __name__ == "__main__":
    rt_args = args.get_cmdarguments()
    try:
        configure_logs(rt_args.loglevel, rt_args.logfile)
        if not rt_args.no_update:
            setup.updates(rt_args.version, rt_args.workers)
    except ValueError as error:
        raise error
    logger.debug("Input args: " + repr(rt_args))
    if rt_args.input is not None:
        logger.info("Input is set, starting calculations")
        out = _uninteractive_session(rt_args)
    else:
        logger.info("No input set, starting interactive session")
        out = _interactive_session(rt_args)
