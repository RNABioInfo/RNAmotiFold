import src.bgap_rna as bgap
import src.args as args
from argparse import Namespace
import src.results
from src.update_motifs import update_hexdumbs
import setup
import configparser
import logging
from pathlib import Path
from typing import Generator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import sys
from typing import Optional


####Logging setup####
import logging

logger = logging.getLogger("RNAmotiFold")


# Interactive session to run multiple predictions in an "interactive" environment
def _interactive_session(
    runtime_arguments: Namespace | configparser.ConfigParser,
) -> list[src.results.algorithm_output | src.results.error | None]:
    """Function is an infinite while Loop that always does one prediction, returns the results and waits for a new input."""
    result_list = []
    if isinstance(runtime_arguments, Namespace):
        proc_obj = bgap.bgap_rna.from_argparse(runtime_arguments)
        output_file = runtime_arguments.output
        pool_boys = runtime_arguments.workers
        csv_separator = runtime_arguments.separator
        name = runtime_arguments.id
    elif isinstance(runtime_arguments, configparser.ConfigParser):
        proc_obj = bgap.bgap_rna.from_config(runtime_arguments, "VARIABLES")
        output_file = runtime_arguments.get("VARIABLES", "output")
        pool_boys = runtime_arguments.getint("VARIABLES", "workers")
        csv_separator = runtime_arguments.get("VARIABLES", "separator")
        name = runtime_arguments.get("VARIABLES", "name")
    else:
        raise TypeError(
            "runtime args is neither a argparse Namespace nor a configparser.Configparser instance. Exiting..."
        )
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
                realtime_input = _input_check(user_input, name)
            except ValueError as error:
                print(error)
            else:
                result = proc_obj.auto_run(
                    realtime_input,
                    o_file=output_file,
                    pool_workers=pool_boys,
                    output_csv_separator=csv_separator,
                )
                result_list.append(result)
    return result_list  # Added result outputting just in case I wanna do something with that down the line.


# Uninteractive session in case of preset input, just does the calculation and exits
def _uninteractive_session(
    runtime_arguments: Namespace | configparser.ConfigParser,
) -> list[src.results.algorithm_output | src.results.error]:
    if isinstance(runtime_arguments, Namespace):
        runtime_input = _input_check(
            runtime_arguments.input,
            runtime_arguments.id,
        )
        proc_obj = bgap.bgap_rna.from_argparse(runtime_arguments)
        output_file = runtime_arguments.output  # type:str
        pool_boys = runtime_arguments.workers  # type: int
        csv_seaparator = runtime_arguments.separator  # type:str
    elif isinstance(runtime_arguments, configparser.ConfigParser):
        runtime_input = _input_check(
            runtime_arguments.get("VARIABLES", "input"),
            runtime_arguments.get("VARIABLES", "name"),
        )
        proc_obj = bgap.bgap_rna.from_config(runtime_arguments)
        output_file = runtime_arguments.get("VARIABLES", "output")  # type:str
        pool_boys = runtime_arguments.getint("VARIABLES", "workers")
        csv_seaparator = runtime_arguments.get("VARIABLES", "separator")  # type:str
    result = proc_obj.auto_run(
        user_input=runtime_input,
        o_file=output_file,
        pool_workers=pool_boys,
        output_csv_separator=csv_seaparator,
    )
    # Added result outputting just in case I wanna do something with that down the line.
    return result


# Finds File type based on file ending
def _find_filetype(file_path: Path) -> None:
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
    SeqIO.FastaIO.FastaIterator
    | SeqIO.QualityIO.FastqPhredIterator
    | Generator[SeqIO.SeqRecord, None, None]
):
    (zipped, filetype) = _find_filetype(file_path)
    if not zipped:
        return SeqIO.parse(file_path, filetype)
    else:
        with gzip.open(file_path, "rt") as handle:
            return SeqIO.parse(handle, filetype)


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
def configure_logs(loglevel: str, logfile: Optional[str]):
    if logfile is not None:
        logging.basicConfig(
            filename=logfile,
            filemode="w+",
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
        if isinstance(rt_args, configparser.ConfigParser):
            configure_logs(
                rt_args.get("VARIABLES", "loglevel"), rt_args.get("VARIABLES", "logfile")
            )
            setup.updates(
                rt_args.getboolean("VARIABLES", "no_update"),
                rt_args.getboolean("VARIABLES", "update_algorithms"),
            )
            user_input = rt_args.get("VARIABLES", "input")
        elif isinstance(rt_args, Namespace):
            configure_logs(rt_args.loglevel, rt_args.logfile)
            setup.updates(rt_args.no_update, rt_args.update_algorithms)
            user_input = rt_args.input
    except ValueError as error:
        raise error

    if user_input is not None:
        logger.info("Input is set, starting calculations")
        out = _uninteractive_session(rt_args)

    else:
        logger.info("No input set, starting interactive session")
        out = _interactive_session(rt_args)
