from __future__ import annotations

import gzip
import logging
from abc import ABC, abstractmethod
from collections import Counter
from collections.abc import Generator
from dataclasses import dataclass
from pathlib import Path
from typing import IO, ClassVar

import Bio.Align
import Bio.SeqIO
import Bio.SeqIO.FastaIO
import Bio.SeqIO.QualityIO
from Bio.SeqRecord import SeqRecord
from RNAmotiFold.bgap_rna.call_handler import MainCall
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


class AlgorithmInput(ABC):
    test_string: ClassVar[str] = ""

    def __init__(self, id, sequence: str):
        self.id: str = id
        self.sequence: str = self._check_sequence(id, sequence)

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
    @abstractmethod
    def from_seqrecord(
        cls, record: Bio.SeqIO.SeqRecord
    ) -> AlgorithmInput: ...

    @classmethod
    @abstractmethod
    def from_alignmentrecord(
        cls, record: Bio.Align.MultipleSeqAlignment, id: str
    ) -> AlgorithmInput: ...

    @classmethod
    @abstractmethod
    def from_strings(
        cls, input_str: str, input_id: str
    ) -> AlgorithmInput: ...


class SingleAlgorithmInput(AlgorithmInput):
    test_string: ClassVar[str] = "AUCGuagcNn"

    def __init__(self, sequence: str, id: str = "N/A"):
        super().__init__(id, sequence)
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
            yield SingleAlgorithmInput.from_seqrecord(obj)

    @classmethod
    def from_alignmentrecord(
        cls, record: Bio.Align.MultipleSeqAlignment, id: str
    ) -> SingleAlgorithmInput:
        raise InputException(
            id,
            "Couldnt create Single Sequence Algorithm input objects from Multiple Sequence alignment",
        )

    @classmethod
    def from_seqrecord(
        cls, record: Bio.SeqIO.SeqRecord
    ) -> SingleAlgorithmInput:
        if record.id is None:
            id = "N/A"
        else:
            id = record.id
        return cls(str(record.seq), id)

    @classmethod
    def from_strings(
        cls, input_str: str, input_id: str
    ) -> SingleAlgorithmInput:
        return cls(sequence=input_str, id=input_id)


class AlignmentAlgorithmInput(AlgorithmInput):
    test_string: ClassVar[str] = "AUCGuacgNn#-_"

    def __init__(self, sequence: str, id: str = "N/A"):
        super().__init__(id, sequence)
        if AlignmentAlgorithmInput.test_equal_lengths(sequence):
            self.input: list[str] = [
                "-f",
                str(TempFileHandler.write_temp_file(self.sequence)),
            ]
        else:
            raise InputException(
                id, "Not all input sequences are of equal length."
            )

    @staticmethod
    def test_equal_lengths(test_string: str, divider: str = "#") -> bool:
        strings = test_string.split(divider)
        return all(len(strings[0]) == len(x) for x in strings)

    @classmethod
    def from_alignmentrecord(
        cls, record: Bio.Align.MultipleSeqAlignment, id: str
    ) -> AlignmentAlgorithmInput:
        CONCATER = "#"
        return cls(CONCATER.join([str(x) for x in record]), id)

    @classmethod
    def _parse_file(
        cls, file_obj: IO, filetype: str
    ) -> Generator[AlignmentAlgorithmInput]:
        alignments = Bio.Align.parse(file_obj, filetype)
        for alignment in alignments:
            yield AlignmentAlgorithmInput.from_alignmentrecord(
                alignment,
                file_obj.name.split("/")[-1].split(".")[
                    0
                ],  # get file name from path
            )

    @classmethod
    def from_seqrecord(
        cls, record: Bio.SeqIO.SeqRecord
    ) -> AlignmentAlgorithmInput:
        if record.id is None:
            id = "N/A"
        else:
            id = record.id
        return cls(str(record.seq), id)

    @classmethod
    def from_strings(
        cls, input_str: str, input_id: str
    ) -> AlignmentAlgorithmInput:
        return cls(sequence=input_str, id=input_id)


@dataclass
class ReadyAlgorithmInput:
    input_obj: SingleAlgorithmInput | AlignmentAlgorithmInput
    base_call: MainCall
    motif_calls: list[MotifCall]

    @property
    def id(self) -> str:
        return self.input_obj.id

    @property
    def call_list(self) -> list[list[str]]:
        if len(self.motif_calls) == 0:
            return [self.base_call.to_list() + self.input_obj.input]

        calls = [
            [*self.base_call.to_list(), *x.call] for x in self.motif_calls
        ]
        full_calls = [x + self.input_obj.input for x in calls]
        return full_calls

    @property
    def call_count(self) -> int:
        if len(self.motif_calls) == 0:
            return 1
        else:
            return len(self.motif_calls)

    @classmethod
    def create_ready_input(
        cls,
        pre_algorithm_input: SingleAlgorithmInput | AlignmentAlgorithmInput,
        base_call: MainCall,
        motif_calls: list[MotifCall],
    ) -> ReadyAlgorithmInput:
        return cls(pre_algorithm_input, base_call, motif_calls)


def get_input(
    base_call: MainCall,
    motif_subcalls: list[MotifCall],
    process_type: ProcessType,
    input: str
    | Bio.SeqIO.SeqRecord
    | Bio.Align.MultipleSeqAlignment
    | Path,
    id: str = "N/A",
):
    inputs = _read_input(base_call, motif_subcalls, process_type, input, id)
    IDs = Counter([x.id for x in inputs])
    duplicate_ids = [key for key in IDs if IDs[key] > 1]
    if any(duplicate_ids):
        raise InputException(
            ", ".join(duplicate_ids),
            "Duplicate IDs found in input, can't process",
        )
    return inputs


def _read_input(
    base_call: MainCall,
    motif_subcalls: list[MotifCall],
    process_type: ProcessType,
    input: str
    | Bio.SeqIO.SeqRecord
    | Bio.Align.MultipleSeqAlignment
    | Path,
    id: str = "N/A",
) -> list[ReadyAlgorithmInput]:
    match process_type:
        case ProcessType.SINGLE:
            AlgorithmInputType = SingleAlgorithmInput
        case ProcessType.ALIGNMENT:
            AlgorithmInputType = AlignmentAlgorithmInput
    match input:
        case str():
            pathd = Path(input)
            if pathd.exists():
                return [
                    ReadyAlgorithmInput.create_ready_input(
                        x, base_call, motif_subcalls
                    )
                    for x in _parse_path(pathd, AlgorithmInputType)
                ]
            else:
                try:
                    return [
                        ReadyAlgorithmInput.create_ready_input(
                            AlgorithmInputType.from_strings(input, id),
                            base_call,
                            motif_subcalls,
                        )
                    ]
                except InputException as exc:
                    raise InputException(
                        id,
                        f"Input is neither a valid RNA sequence nor an existing path: {input!r}",
                    ) from exc
        case Bio.SeqIO.SeqRecord():
            return [
                ReadyAlgorithmInput.create_ready_input(
                    AlgorithmInputType.from_seqrecord(input),
                    base_call,
                    motif_subcalls,
                )
            ]
        case Bio.Align.MultipleSeqAlignment():
            return [
                ReadyAlgorithmInput.create_ready_input(
                    AlgorithmInputType.from_alignmentrecord(input, id),
                    base_call,
                    motif_subcalls,
                )
            ]
        case Path():
            return [
                ReadyAlgorithmInput.create_ready_input(
                    x, base_call, motif_subcalls
                )
                for x in _parse_path(input, AlgorithmInputType)
            ]


def _parse_path(
    user_input: Path,
    AlgorithmInputType: type[
        SingleAlgorithmInput | AlignmentAlgorithmInput
    ],
) -> list[SingleAlgorithmInput | AlignmentAlgorithmInput]:
    if user_input.is_file():
        return _read_file(user_input, AlgorithmInputType)
    elif user_input.is_dir():
        return _read_dir(user_input, AlgorithmInputType)
    else:
        raise OSError(
            "Given Path could not be identified as either a directory or a file"
        )


def _read_file(
    user_input: Path,
    AlgorithmInputType: type[
        SingleAlgorithmInput | AlignmentAlgorithmInput
    ],
) -> list[SingleAlgorithmInput | AlignmentAlgorithmInput]:
    (zipped, potential_filetypes) = _find_filetype(user_input)
    if zipped:
        open_funct = gzip.open
    else:
        open_funct = open
    with open_funct(user_input, "rt") as handle:
        for filetype in potential_filetypes:
            handle.seek(0)
            try:
                inputs = list(
                    AlgorithmInputType._parse_file(handle, filetype)
                )
            except (TypeError, OSError, ValueError):
                logger.info(
                    f"Filetype {filetype} couldnt be used to parse file {user_input!s}, trying other options..."
                )
                continue
            if inputs:
                return inputs
    raise InputException(
        str(user_input),
        f"Couldn't be processed using recognized filetypes: {', '.join(potential_filetypes)}",
    )


def _read_dir(
    user_input: Path,
    AlgorithmInputType: type[
        SingleAlgorithmInput | AlignmentAlgorithmInput
    ],
) -> list[SingleAlgorithmInput | AlignmentAlgorithmInput]:
    inputs: list = []
    for child in user_input.rglob("*.*"):
        if child.is_file():
            try:
                inputs.extend(_read_file(child, AlgorithmInputType))
            except InputException as e:
                logger.warning(e)
                continue
    return inputs


def _find_filetype(file_path: Path) -> tuple[bool, list[str]]:
    """Function for determining file type and zip status of user input, if file suffix is unknown this will try to use the suffix as a argument for BioSeqIO parse"""
    try:
        if not file_path.exists():
            raise FileNotFoundError(
                f"Could not find set input file {file_path}"
            )
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
