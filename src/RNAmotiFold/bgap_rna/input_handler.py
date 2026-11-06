from __future__ import annotations

import gzip
import logging
from abc import ABC, abstractmethod
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from typing import IO, ClassVar, Self

import Bio.Align
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO
from Bio.SeqRecord import SeqRecord
from RNAmotiFold.bgap_rna.motif_handler import MotifCall
from RNAmotiFold.bgap_rna.tempfile_handler import TempFileHandler
from RNAmotiFold.Misc.EnumTypes import ProcessType

logger = logging.getLogger(__name__)


class InputException(Exception):
    def __init__(
        self,
        input_obj_id: str,
        msg="Exception raised during creation of Algorithm Inputs",
    ):
        self.input_obj = input_obj_id
        self.msg = msg
        super().__init__(self.msg)

    def __str__(self):
        return f"Error during processing of {self.input_obj}: {self.msg}"


def create_input_handler(process_type: ProcessType) -> InputHandler:
    match process_type:
        case ProcessType.ALIGNMENT:
            return InputHandler(process_type, AlignmentAlgorithmInput)
        case ProcessType.SINGLE:
            return InputHandler(process_type, SingleAlgorithmInput)


class InputHandler[T: (AlgorithmInput)]:
    def __init__(self, process_type: ProcessType, input_type: type[T]):
        self.process_type = process_type
        self.input_type = input_type

    def _read_file(self, user_input: Path) -> list[T]:
        (zipped, potential_filetypes) = InputHandler._find_filetype(user_input)
        if zipped:
            open_funct = gzip.open
        else:
            open_funct = open
        with open_funct(user_input, "rt") as handle:
            for filetype in potential_filetypes:
                handle.seek(0)
                try:
                    return list(self.input_type._parse_file(handle, filetype))
                except (TypeError, OSError):
                    logger.info(
                        f"Filetype {filetype} couldnt be used to parse input file, trying other options..."
                    )
                    continue
        raise InputException(
            f"{user_input!s} couldn't be processed. Recognized potential filetypes as {', '.join(potential_filetypes)}, but couldnt create input objects"
        )

    def _read_dir(self, user_input: Path) -> list[T]:
        inputs: list[T] = []
        for child in user_input.iterdir():
            if child.is_file():
                inputs.extend(self._read_file(child))
        return inputs

    def _parse_path(self, user_input: Path) -> list[T]:
        if user_input.is_file():
            return self._read_file(user_input)
        elif user_input.is_dir():
            return self._read_dir(user_input)
        else:
            raise OSError(
                "Given Path could not be identified as either a directory or a file"
            )

    def _read_str_input(self, input: str, id: str = "N/A") -> list[T]:
        return [self.input_type.from_strings(input, id)]

    def _read_seqrecord(self, input_record: Bio.SeqIO.SeqRecord):
        return self.input_type.from_seqrecord(input_record)

    def _read_msa(self, input_record: Bio.Align.MultipleSeqAlignment, id: str = "N/A"):
        return self.input_type.from_alignmentrecord(input_record, id)

    def read_input(
        self,
        input: str | Bio.SeqIO.SeqRecord | Bio.Align.MultipleSeqAlignment | Path,
        id: str = "N/A",
    ) -> list[T]:
        match input:
            case str():
                pathd = Path(input)
                if pathd.exists():
                    return self._parse_path(pathd)
                else:
                    try:
                        return self._read_str_input(input, id)
                    except InputException as exc:
                        raise InputException(
                            id,
                            f"Input is neither a valid RNA sequence nor an existing path: {input!r}",
                        ) from exc
            case Bio.SeqIO.SeqRecord():
                return [self._read_seqrecord(input)]
            case Bio.Align.MultipleSeqAlignment():
                return [self._read_msa(input)]
            case Path():
                return self._parse_path(input)

    # Finds File type based on file ending
    @staticmethod
    def _find_filetype(file_path: Path) -> tuple[bool, list[str]]:
        """Function for determining file type and zip status of user input, if file suffix is unknown this will try to use the suffix as a argument for BioSeqIO parse"""
        try:
            if not file_path.exists():
                raise FileNotFoundError(f"Could not find set input file {file_path}")
            if file_path.suffixes[-1] == ".gz":
                file_extension = file_path.suffixes[-2]
                input_zipped = True
            else:
                file_extension = file_path.suffixes[-1]
                input_zipped = False
        except IndexError:
            raise InputException(
                str(file_path),
                "Couldn't correctly access file endings, please add appropriate suffix for you filetype",
            )
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


class AlgorithmInput(ABC):
    test_string: ClassVar[str] = ""

    def __init__(self, id, sequence):
        self.id: str = id
        self.sequence: str = sequence
        self.input: list[str] = [sequence]

    @classmethod
    def from_strings(cls, input_str: str, input_id: str) -> Self:
        return cls(sequence=input_str, id=input_id)

    def _check_sequence(self, id: str, sequence: str) -> str:
        sequence = sequence.replace("t", "u")
        sequence = sequence.replace("T", "U")
        sequence = sequence.replace("-", "_")
        test_list = [
            (idx, char)
            for idx, char in enumerate(sequence)
            if char not in self.test_string
        ]
        if len(test_list) > 0:
            raise InputException(
                id,
                f"Found irregular characters in sequence {id}: {','.join([f'Position {elem[0]}:{elem[1]}' for elem in test_list])}",
            )
        return sequence

    @classmethod
    @abstractmethod
    def _parse_file(cls, file_obj: IO, filetype: str) -> Generator: ...

    @classmethod
    def from_seqrecord(cls, record: Bio.SeqIO.SeqRecord):
        if record.id is None:
            id = "N/A"
        else:
            id = record.id
        return cls(str(record.seq), id)

    @classmethod
    @abstractmethod
    def from_alignmentrecord(
        cls, record: Bio.Align.MultipleSeqAlignment, id: str
    ) -> Self: ...


class SingleAlgorithmInput(AlgorithmInput):
    test_string: ClassVar[str] = "AUCGuagcNn"

    def __init__(self, sequence: str, id: str = "N/A"):
        self.id = id
        self.sequence = self._check_sequence(id, sequence)
        self.input: list[str] = [self.sequence]

    @classmethod
    def _parse_file(
        cls, file_obj: IO, filetype: str
    ) -> Generator[SingleAlgorithmInput]:
        objs: (
            Bio.SeqIO.FastaIO.FastaIterator
            | Bio.SeqIO.QualityIO.FastqPhredIterator
            | Generator[SeqRecord, None, None]
        ) = Bio.SeqIO.parse(file_obj, filetype)
        for obj in objs:
            try:
                yield SingleAlgorithmInput.from_seqrecord(obj)
            except InputException as e:
                logger.error(e)

    @classmethod
    def from_alignmentrecord(
        cls, record: Bio.Align.MultipleSeqAlignment, id: str
    ) -> Self:
        raise InputException(
            id,
            "Couldnt create Algorithm input objects from Multiple Sequence alignment",
        )


class AlignmentAlgorithmInput(AlgorithmInput):
    test_string: ClassVar[str] = "AUCGuacgNn#-_"

    def __init__(self, sequence: str, id: str = "N/A"):
        self.id = id
        self.sequence = self._check_sequence(id, sequence)
        if AlignmentAlgorithmInput.test_equal_lengths(sequence):
            self.input: list[str] = [
                "-f",
                str(TempFileHandler.write_temp_file(self.sequence)),
            ]
        else:
            raise InputException(id, "Not all input sequences are of equal length.")

    @staticmethod
    def test_equal_lengths(test_string: str, divider: str = "#") -> bool:
        strings = test_string.split(divider)
        return all(len(strings[0]) == len(x) for x in strings)

    @classmethod
    def from_alignmentrecord(
        cls, record: Bio.Align.MultipleSeqAlignment, id: str
    ) -> Self:
        CONCATER = "#"
        return cls(CONCATER.join(record), id)

    @classmethod
    def _parse_file(
        cls, file_obj: IO, filetype: str
    ) -> Generator[AlignmentAlgorithmInput]:
        alignments = Bio.Align.parse(file_obj, filetype)
        for alignment in alignments:
            try:
                yield AlignmentAlgorithmInput.from_alignmentrecord(
                    alignment,
                    file_obj.name.split("/")[-1].split(".")[
                        0
                    ],  # get file name from path
                )
            except InputException as e:
                logger.warning(e)


@dataclass
class ReadyAlgorithmInput[T: AlgorithmInput]:
    input_obj: T
    base_call: list[str]
    motif_calls: list[MotifCall]

    @property
    def id(self) -> str:
        return self.input_obj.id

    @property
    def call_list(self) -> list[list[str]]:
        if len(self.motif_calls) == 0:
            return [self.base_call + self.input_obj.input]
        calls = []
        for motif_call in self.motif_calls:
            calls.append(self.base_call + motif_call.call + self.input_obj.input)
        return calls

    @property
    def call_count(self) -> int:
        if len(self.motif_calls) == 0:
            return 1
        else:
            return len(self.motif_calls)

    @classmethod
    def create_ready_input(
        cls, pre_algorithm_input: T, base_call: list[str], motif_calls: list[MotifCall]
    ) -> ReadyAlgorithmInput[T]:
        return cls(pre_algorithm_input, base_call, motif_calls)
