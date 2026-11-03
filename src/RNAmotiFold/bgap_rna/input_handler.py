import glob
import gzip
import logging
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Self

import Bio.Align
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO
from Bio.SeqRecord import SeqRecord
from RNAmotiFold.bgap_rna.tempfile_handler import TempFileHandler
from RNAmotiFold.Misc.EnumTypes import ProcessType

logger = logging.getLogger(__name__)


# List flattening
def flatten[T](xss: list[list[T]]) -> list[T]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]


@dataclass
class AlgorithmInput:
    _id: str | None
    input_str: str
    process_type: ProcessType
    _call: str = ""

    @property
    def runtime_call(self) -> str:
        """Creates runtime call of object by combining its call and input_str"""
        if self.call == "":
            raise ValueError("No call set yet for this algorithm input object")
        if self.process_type == ProcessType.ALIGNMENT:
            tmp = TempFileHandler().write_temp_file(self.input_str)
            return f"{self.call} -f {tmp.name}"
        else:
            return f"{self.call} {self.input_str}"

    @property
    def call(self) -> str:
        return self._call

    @call.setter
    def call(self, new_call: str) -> None:
        self._call = new_call

    @property
    def id(self) -> str:
        if self._id is None:
            return "N/A"
        else:
            return self._id

    @classmethod
    def from_generators(
        cls,
        process_type: ProcessType,
        input_generator: Any,
    ) -> list["AlgorithmInput"]:
        """Classmethod to create list of algroithm inputs (concated alignment sequences or individual mfe prediction sequences), depending on the process type
        given to the input handler by the bgap object.
        """
        returnlist: list[AlgorithmInput] = []
        match process_type:
            case ProcessType.ALIGNMENT:
                for alignment in input_generator:
                    for seq in alignment:
                        if any(c not in "NAUCGTnaucgt+_-#" for c in set(seq)):
                            raise ValueError(
                                "Input was neither a viable file path nor a viable RNA or DNA sequence"
                            )
                    if all(len(item) == len(alignment[0]) for item in alignment):
                        concat: str = "#".join(alignment)
                        returnlist.append(
                            AlgorithmInput(
                                input_generator.source.name.split(".")[0],
                                concat.replace("T", "U")
                                .replace("-", "_")
                                .replace("t", "u"),
                                process_type=process_type,
                            )
                        )
                    else:
                        raise ValueError(
                            "Not all sequences have the same length in your alignment, can't process input"
                        )
            case ProcessType.SINGLE:
                for sequence in input_generator:
                    returnlist.append(
                        AlgorithmInput(
                            _id=sequence.id,
                            input_str=str(sequence.seq),
                            process_type=process_type,
                        )
                    )
        return returnlist


class InputHandler:
    """Class to manage inputs, takes the chosen algorithm, the given inputs etc into account and creates calls from them"""

    def __iter__(self) -> Self:
        return self

    def __next__(self) -> AlgorithmInput:
        if self._index < len(self.user_input):
            item = self.user_input[self._index]
            self._index += 1
            return item
        else:
            self._index = 0
            raise StopIteration

    def __init__(
        self,
        process_type: ProcessType,
        user_input: str | None,
    ) -> None:
        self.process_type: ProcessType = process_type
        if user_input is not None:
            self.user_input: list[AlgorithmInput] = InputHandler.read_input(
                process_type, user_input
            )
        self._index = 0
        logger.debug(f"Created input handler: {vars(self)}")

    @staticmethod
    def read_input(
        process_type: ProcessType,
        user_input: str,
        id: str = "N/A",
    ) -> list[AlgorithmInput]:
        """Read input given by user and tries to determine if it's a file, a directory or just a single sequence.
        I really wanted to strictly type this but the number of different possible iterators make it basically impossible.
        This will always return a list of Iterators, even if there is only on MSA or seq in the file, we take care of that after."
        """
        logger.debug("Attempting to read input")
        pathd = Path(
            user_input.strip()
        )  # If this one doesn't work then it's also not a directory
        if pathd.resolve().is_file():
            logger.debug("Input recognized as file, reading now")
            result = [InputHandler._read_input_file(pathd, process_type)]
            logger.debug("Input file read worked, generating algorithm inputs")
            return flatten(
                [AlgorithmInput.from_generators(process_type, x) for x in result]
            )
        elif pathd.resolve().is_dir():
            logger.debug("Input recognized as folder, parsing now")
            result = [
                InputHandler._read_input_file(pathd / Path(file), process_type)
                for file in glob.glob("*", root_dir=pathd)
            ]
            logger.debug("Folder parsing successfull, generating algorithm inputs")
            return flatten(
                [AlgorithmInput.from_generators(process_type, x) for x in result]
            )
        # Wasn't a file or a directory, so it's probably a sequence
        if process_type == "ali":
            teststring: str = "NAUCGTnaucgt_#"
            if not InputHandler.test_equal_lengths(user_input):
                raise ValueError(
                    "Not all strings in your alignment have the same lengths, can't process"
                )
        else:
            teststring: str = "NAUCGTnaucgt"
        if any(c not in teststring for c in set(user_input.strip())):
            raise ValueError(
                "Input was neither a viable file path nor a viable RNA or DNA sequence"
            )

        # No error was raised so we can assume that it is not a directory or file in the system but it is a viable sequence
        return [
            AlgorithmInput(
                _id=id,
                input_str=user_input.strip(),
                process_type=process_type,
            )
        ]

    @staticmethod
    def test_equal_lengths(test_string: str, divider: str = "#") -> bool:
        strings = test_string.split(divider)
        base_length = len(strings[0])
        return all(len(x) == base_length for x in strings)

    @staticmethod
    def script_read_input(
        process_type: ProcessType,
        user_input: (str | SeqRecord | Path | list[SeqRecord] | list[str]),
        id: str = "N/A",
    ) -> list[AlgorithmInput] | list[Any]:
        match user_input:
            case str():
                return InputHandler.read_input(process_type, user_input, id)
            case Path():
                return InputHandler.read_input(process_type, str(user_input), id)
            case SeqRecord():
                return [
                    AlgorithmInput(
                        _id=user_input.id,
                        input_str=str(user_input.seq),
                        process_type=process_type,
                    )
                ]
            case list():
                if all(isinstance(x, str) for x in user_input):
                    return flatten(
                        [
                            InputHandler.read_input(process_type, x, id)  # type: ignore
                            for x in user_input
                        ]
                    )
                elif all(isinstance(x, SeqRecord) for x in user_input):
                    return [
                        AlgorithmInput(
                            _id=x.id,  # type: ignore
                            input_str=str(x.seq),  # type: ignore
                            process_type=process_type,
                        )
                        for x in user_input
                    ]
                else:
                    raise ValueError("Unknown type in input detected")
            case _:
                raise ValueError("Couldnt parse input")

    # Read input file
    @staticmethod
    def _read_input_file(
        file_path: Path, process_type: ProcessType
    ) -> (
        Bio.SeqIO.FastaIO.FastaIterator
        | Bio.SeqIO.QualityIO.FastqPhredIterator
        | Generator[SeqRecord, None, None]
    ):
        zipped, filetype = InputHandler._find_filetype(file_path)
        if (
            process_type == ProcessType.ALIGNMENT
        ):  # Either parse it with the Align parser if we're doing alignments or with the SeqIO if we're doing MFE/PFC calculations
            parse_func = Bio.Align.parse
        else:
            parse_func = Bio.SeqIO.parse
        if not zipped:
            for option in filetype:
                logger.debug(
                    f"Recognized input {file_path} as not compressed, reading as {option} file."
                )
                return parse_func(file_path, option)

        else:
            with gzip.open(file_path, "rt") as handle:
                for option in filetype:
                    logger.debug(
                        f"Recognized input {file_path} as compressed, decrompressing and reading as {option} file"
                    )
                    return parse_func(handle, option)
        raise NotImplementedError(
            "Could not read given input file, probably because the given file type is not implemented. Try a different one or write up an issue on GitHub!"
        )

    # Finds File type based on file ending
    @staticmethod
    def _find_filetype(file_path: Path) -> tuple[bool, list[str]]:
        """Function for determining file type and zip status of user input, if file suffix is unknown this will try to use the suffix as a argument for BioSeqIO parse"""
        if file_path.suffixes[-1] == ".gz":
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
                filetype = ["fasta", "emboss", "tabular"]
            case ".fastq" | ".fq":
                filetype = ["fastq"]
            case ".stk" | ".stockholm" | ".sto":
                filetype = ["stockholm"]
            case ".phy":
                filetype = ["phylip"]
            case _:
                filetype = [file_extension]
        return (input_zipped, filetype)
