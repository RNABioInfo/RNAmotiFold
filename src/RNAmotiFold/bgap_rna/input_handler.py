from typing import Literal, Any, ClassVar
from dataclasses import dataclass
from _collections_abc import Generator
from pathlib import Path
import Bio.SeqIO.FastaIO, Bio.SeqIO, Bio.SeqIO.QualityIO
import Bio.Align
from Bio.SeqRecord import SeqRecord
import gzip
import glob
import tempfile
from os import remove
import logging

logger = logging.getLogger(__name__)


# List flattening
def flatten(xss: list[list[Any]]) -> list[Any]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]


@dataclass
class algorithm_input:

    tmp_files: ClassVar[list[Path]] = []

    id: str
    input_str: str
    process_type: Literal["ali", "single"]
    _call: str = ""

    @staticmethod
    def cleanup_temps():
        """Removes temporary files created for RNAmotiAlign"""
        logger.debug(f"Cleaning up input temp files {algorithm_input.tmp_files}")
        for file in algorithm_input.tmp_files:
            remove(file)

    @property
    def runtime_call(self) -> str:
        """Creates runtime call of object by combining its call and input_str"""
        if self.call == "":
            raise ValueError("No call set yet for this algorithm input object")
        if self.process_type == "ali":
            logger.info("Writing alignment folding input to temp file")
            tmp = tempfile.NamedTemporaryFile(delete=False, delete_on_close=False)
            tmp.write(self.input_str.encode())
            logger.info(f"Input for id {self.id} written to {tmp.name}")
            algorithm_input.tmp_files.append(Path(tmp.name))
            return " ".join([self.call, f"-f {tmp.name}"])
        else:
            return " ".join([self.call, self.input_str])

    @property
    def call(self):
        return self._call

    @call.setter
    def call(self, new_call: str):
        self._call = new_call

    @classmethod
    def from_generators(
        cls, process_type: Literal["ali", "single"], input_generator: Any
    ) -> list["algorithm_input"]:
        """Classmethod to create list of algroithm inputs (concated alignment sequences or individual mfe prediction sequences), depending on the process type
        given to the input handler by the bgap object.
        """
        returnlist: list[algorithm_input] = []
        match process_type:
            case "ali":
                for alignment in input_generator:
                    try:
                        seqs = [
                            x.replace("-", "_").replace("T", "U").replace("t", "u").upper()
                            for x in alignment
                        ]
                        for seq in seqs:
                            if any(c not in "NAUCGTnaucgt+_-#" for c in set(seq)):
                                raise ValueError(
                                    "Input was neither a viable file path nor a viable RNA or DNA sequence"
                                )
                    except:
                        raise TypeError("Couldn't convert sequences from DNA to RNA")
                    concat = "#".join(seqs)
                    returnlist.append(algorithm_input("N/A", concat, process_type=process_type))
            case "single":
                for sequence in input_generator:
                    returnlist.append(
                        algorithm_input(
                            id=sequence.id, input_str=str(sequence.seq), process_type=process_type
                        )
                    )
            case _:
                raise ValueError("Unknown state detected during algorithm input generation")
        return returnlist


class input_handler:
    """Class to manage inputs, takes the chosen algorithm, the given inputs etc into account and creates calls from them"""

    def __iter__(self):
        return self

    def __next__(self):
        if self._index < len(self.user_input):
            item = self.user_input[self._index]
            self._index += 1
            return item
        else:
            self._index = 0
            raise StopIteration

    def __init__(self, process_type: Literal["single", "ali"], user_input: str | None) -> None:
        self.process_type: Literal["single"] | Literal["ali"] = process_type
        if user_input is not None:
            self.user_input: list[algorithm_input] = input_handler.read_input(
                process_type, user_input
            )
        self._index = 0
        logger.debug(f"Created input handler: {vars(self)}")

    @staticmethod
    def read_input(
        process_type: Literal["single", "ali"], user_input: str, id: str = "N/A"
    ) -> list[algorithm_input]:
        """Read input given by user and tries to determine if it's a file, a directory or just a single sequence.
        I really wanted to strictly type this but the number of different possible iterators make it basically impossible.
        This will always return a list of Iterators, even if there is only on MSA or seq in the file, we take care of that after."
        """
        logger.debug("Attempting to read input")
        try:
            pathd = Path(
                user_input.strip()
            )  # If this one doesn't work then it's also not a directory
            if pathd.resolve().is_file():
                logger.debug("Input recognized as file, reading now")
                result = [input_handler._read_input_file(pathd, process_type)]
                logger.debug("Input file read worked, generating algorithm inputs")
                return flatten([algorithm_input.from_generators(process_type, x) for x in result])
            elif pathd.resolve().is_dir():
                logger.debug("Input recognized as folder, parsing now")
                result = [
                    input_handler._read_input_file(pathd / Path(file), process_type)
                    for file in glob.glob("*", root_dir=pathd)
                ]
                logger.debug("Folder parsing successfull, generating algorithm inputs")
                return flatten([algorithm_input.from_generators(process_type, x) for x in result])
            if any(c not in "NAUCGTnaucgt+_-#" for c in set(user_input.strip())):
                raise ValueError(
                    "Input was neither a viable file path nor a viable RNA or DNA sequence"
                )
        except OSError as e:
            raise e
        except ValueError as e:
            raise e
        # No error was raised so we can assume that it is not a directory or file in the system but it is a viable sequence
        return [algorithm_input(id=id, input_str=user_input.strip(), process_type=process_type)]

    # Read input file
    @staticmethod
    def _read_input_file(
        file_path: Path, process_type: Literal["single", "ali"]
    ) -> (
        Bio.SeqIO.FastaIO.FastaIterator
        | Bio.SeqIO.QualityIO.FastqPhredIterator
        | Generator[SeqRecord, None, None]
    ):
        zipped, filetype = input_handler._find_filetype(file_path)
        if (
            process_type == "ali"
        ):  # Either parse it with the Align parser if we're doing alignments or with the SeqIO if we're doing MFE/PFC calculations
            parse_func = Bio.Align.parse  # type: ignore
        else:
            parse_func = Bio.SeqIO.parse  # type: ignore
        if not zipped:
            logger.debug(f"Recognized input {file_path} as not compressed, reading as {filetype} file.")
            for option in filetype:
                try:
                    return parse_func(file_path, option)  # type: ignore Both of these ignores are because of funky typing on parse from Bio
                except:
                    pass
        else:
            logger.debug(f"Recognized input {file_path} as compressed, decrompressing and reading as {filetype} file")
            with gzip.open(file_path, "rt") as handle:
                for option in filetype:
                    try:
                        return parse_func(handle, option)  # type: ignore
                    except:
                        pass
        raise NotImplementedError(
            "Could not read given input file, probably because the given file type is not implemented. Try a different one or write up an issue on GitHub!"
        )

    # Finds File type based on file ending
    @staticmethod
    def _find_filetype(file_path: Path) -> tuple[bool, list[str]]:
        """Function for determining file type and zip status of user input, if file suffix is unknown this will try to use the suffix as a argument for BioSeqIO parse"""
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
                filetype = ["fasta", "emboss", "tabular"]
            case ".fastq" | ".fq":
                filetype = ["fastq"]
            case ".stk" | ".stockholm" | ".sto":
                filetype = ["stockholm"]
            case "phy":
                filetype = ["phylip"]
            case _:
                filetype = [file_extension]
        return (input_zipped, filetype)