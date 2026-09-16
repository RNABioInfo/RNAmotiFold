from contextlib import redirect_stdout
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import FastaIO
import tempfile
import subprocess
from os import path, remove
from sys import stdout
from ..results import algorithm_output, error
from typing import Optional, Any
from pathlib import Path


class AlignmentProcess:
    """Class to handle alignment folding process, including input formatting and running the bgap-algorithm as a subprocess."""

    concater: str = "#"

    def __init__(self):
        pass

    @staticmethod
    def check_sequence(sequence: str):
        """Checks if the input sequence is valid, returns True if valid, raises ValueError if not"""
        if len(sequence) == 0:
            raise ValueError("Input sequence is empty.")
        valid_nucleotides = set("ACGUacguNn-")
        if not all(nucleotide in valid_nucleotides for nucleotide in sequence):
            raise ValueError(
                f"Invalid characters found in the input sequence: {sequence}. Only A, C, G, U, N and - are allowed."
            )
        return True

    @staticmethod
    def format_input(
        input_sequences: list[SeqRecord] | FastaIO.FastaIterator,
    ):
        """Formats input to a bgap-algorithm readable string and returns it. Concats with given concater, 
        replaces gap - with _, replaces T with U, checks input as well for unknown characters
        """
        if not isinstance(input_sequences, list):
            input_sequences = list(input_sequences)
        try:
            formatted_sequences = [
                str(x.seq).replace("-", "_").replace("T", "U").upper()
                for x in input_sequences
                if AlignmentProcess.check_sequence(str(x.seq))
            ]
        except ValueError as e:
            raise ValueError(
                "Input sequences could not be formatted to a bgap-algorithm readable string. Error: "
                + str(e)
            )
        return AlignmentProcess.concater.join(formatted_sequences)

    @staticmethod
    def _write_formatted_seq_to_temp(formatted_input: str):
        """Writes the formatted input to a temporary file and returns the path to that file"""
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".tmp", prefix="bgap_input_")
        with open(temp_file.name, "w") as f:
            f.write(formatted_input)
        return temp_file.name
    
    @staticmethod
    def run(
        run_id: str,
        call: str,
        input_sequences: list[SeqRecord] | FastaIO.FastaIterator,
        output_file: Optional[str | Path] = None,
    ) -> list[algorithm_output | Any]:
        """Runs the bgap-algorithm with the given call and input, returns the output as a string"""
        # Format Sequence inputs
        ready_input = AlignmentProcess.format_input(input_sequences)
        # Write them to a temp file to use as input for alignment folding
        temp_input_path = AlignmentProcess._write_formatted_seq_to_temp(ready_input)
        # run the the alignment folding algorithm as a subprocess
        try:
            subproc_out = subprocess.run(
                [call, f"-f {temp_input_path}"], capture_output=True, text=True, check=True
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Alignment folding algorithm failed with error: {e.stderr}")
        else:
            if not subproc_out.returncode:
                return_value = algorithm_output(
                    name=run_id, result_str=subproc_out.stdout, stderr=[subproc_out.stderr]
                )
                if output_file is not None:
                    with open(output_file, "a+") as write_file:
                        with redirect_stdout(write_file):
                            return_value.write_results(initiated=False)
                else:
                    return_value.write_results(initiated=False)
                    stdout.flush()
            else:
                return_value = error(id=run_id, error=subproc_out.stderr)
        finally:
            # Clean up the temporary input file using OS functions
            if path.exists(temp_input_path):
                remove(temp_input_path)
        return [return_value]
