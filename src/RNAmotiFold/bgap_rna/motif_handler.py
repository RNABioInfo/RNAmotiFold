from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Self

from RNAmotiFold.bgap_rna.input_handler import flatten
from RNAmotiFold.bgap_rna.tempfile_handler import TempFileHandler
from RNAmotiFold.input.parameters import ScriptParameters
from RNAmotiFold.Misc.EnumTypes import LoopType, MotifOrigin

from RNAmotiFold import RNALOOPS_PATH, RNAMOTIFOLD_ROOT_DIR

logger = logging.getLogger(__name__)


def get_motif_file_paths(version: str) -> list[Path]:
    return [
        x
        for x in RNALOOPS_PATH.joinpath(
            "Misc",
            "Applications",
            "RNAmotiFold",
            "motifs",
            "versions",
            "separated",
            f"{version}_separated",
        ).rglob("*.csv")
    ]


@dataclass
class MotifSet:
    motif_list: list[Motif]
    abbreviation: str
    loop_type: LoopType

    @classmethod
    def from_list(cls, motlist: list[Motif]):
        types = {x.loop_type for x in motlist}
        if len(types) > 1:
            raise ValueError("Motifs of different loop types grouped together")
        abbreviations = {x.abbreviation for x in motlist}
        if len(abbreviations) > 1:
            raise ValueError("Motifs of different abberviations grouped together")
        return cls(
            motif_list=motlist, abbreviation=abbreviations.pop(), loop_type=types.pop()
        )

    def write_self(self) -> None:
        if len(self.motif_list) > 0:
            self.filepath = TempFileHandler.write_temp_file(
                "\n".join([x.csv_format for x in self.motif_list])
            )
        else:
            raise ValueError(
                f"Didn't write file for {self.abbreviation}, motif list empty."
            )


@dataclass
class Motif:
    """Class for storing the information about a motif, including its abbreviation, the sequence and the type of motif it is"""

    sequence: str
    abbreviation: str
    loop_type: LoopType
    origin: MotifOrigin

    @property
    def csv_format(self) -> str:
        return self.sequence + "," + self.abbreviation

    @classmethod
    def from_string(
        cls, motif_string: str, loop_type: LoopType, origin: MotifOrigin
    ) -> Self:
        """Creates a Motif object from a string, the string should be in the format of 'sequence,abbreviation'"""
        split = motif_string.split(",")
        if len(split) != 2:
            raise ValueError(
                f"Motif string {motif_string} is not in the correct format, should be 'sequence,abbreviation'"
            )
        return cls(split[0].strip(), split[1].strip(), loop_type, origin)

    @classmethod
    def from_file(
        cls, file_path: Path, loop_type: LoopType, origin: MotifOrigin
    ) -> Iterator[Self]:
        """Creates a Motif object from a file, the file should be in the format of 'sequence,abbreviation' per line"""
        with open(file_path, "r") as motif_file:
            for line in motif_file:
                yield cls.from_string(line.strip(), loop_type, origin)


@dataclass
class MotifCall:
    hairpin_path: Path | None
    internal_path: Path | None
    bulge_path: Path | None

    @property
    def call(self):
        returnlist: list[str] = []
        if self.hairpin_path is not None:
            returnlist.append(f"-X {self.hairpin_path!s}")
            returnlist.append("-L 1")
        if self.internal_path is not None:
            returnlist.append(f"-Y {self.internal_path!s}")
            returnlist.append("-E 1")
        if self.bulge_path is not None:
            returnlist.append(f"-Z {self.bulge_path!s}")
            returnlist.append("-G 1")
        return returnlist


class MotifHandler:
    """Motif Handler Class, creating MotifCall objects, each of which represents a unique Motif Call."""

    def __init__(
        self,
        motifs: str | None = None,
        single_motif_mode: bool = False,
        custom_hairpin_filepath: Path | None = None,
        custom_internal_filepath: Path | None = None,
        custom_bulge_filepath: Path | None = None,
        replace_hairpins: bool = False,
        replace_internals: bool = False,
        replace_bulges: bool = False,
        version: str = "4_14",
    ) -> None:
        custom_hairpins = MotifHandler.get_custom_motifs(
            custom_hairpin_filepath, LoopType.HAIRPIN, MotifOrigin.CUSTOM
        )
        custom_internals = MotifHandler.get_custom_motifs(
            custom_internal_filepath, LoopType.INTERNAL, MotifOrigin.CUSTOM
        )
        custom_bulges = MotifHandler.get_custom_motifs(
            custom_bulge_filepath, LoopType.BULGE, MotifOrigin.CUSTOM
        )
        default_motifs = flatten(
            [
                list(
                    Motif.from_file(
                        x,
                        LoopType(str(x.parent.name[:-1])),
                        origin=MotifOrigin.DEFAULT,
                    )
                )
                for x in get_motif_file_paths(version)
            ]
        )
        motifs_list: list[Motif] = (
            default_motifs + custom_hairpins + custom_internals + custom_bulges
        )
        if replace_hairpins:
            motifs_list = [
                x
                for x in motifs_list
                if not (
                    x.loop_type == LoopType.HAIRPIN and x.origin == MotifOrigin.DEFAULT
                )
            ]
        if replace_internals:
            motifs_list = [
                x
                for x in motifs_list
                if not (
                    x.loop_type == LoopType.INTERNAL and x.origin == MotifOrigin.DEFAULT
                )
            ]
        if replace_bulges:
            motifs_list = [
                x
                for x in motifs_list
                if not (
                    x.loop_type == LoopType.BULGE and x.origin == MotifOrigin.DEFAULT
                )
            ]
        if not motifs: #FIXME something doesnt quite work with the motif strings rn fix tomorrow, this is a temporary solution
            self.motifs = motifs_list
        else:
            self.motifs = [x for x in motifs_list if x.abbreviation in motifs]
        self.single_motif_mode = single_motif_mode

    @staticmethod
    def get_custom_motifs(
        custom_path: Path | None, motif_type: LoopType, origin: MotifOrigin
    ) -> list[Motif]:
        """Returns a list of Motif objects from a given file path"""
        if custom_path is None:
            return []
        return list(Motif.from_file(custom_path, motif_type, origin))

    @classmethod
    def from_script_parameters(cls, params: ScriptParameters) -> Self:
        return cls(
            params.motif_list,
            params.single_motif_mode,
            params.custom_hairpins,
            params.custom_internals,
            params.custom_bulges,
            params.replace_hairpins,
            params.replace_internals,
            params.replace_bulges,
            params.version,
        )

    def create_motif_calls(self) -> list[MotifCall]:
        # New implementation pre-filters motifs, so self.motifs only contains the motifs we want.
        # First case: Combinatorial where we can just put everything with the same motif type together
        if not self.single_motif_mode:
            if all(x.origin == MotifOrigin.DEFAULT for x in self.motifs):
                return []  # Default behavior, no custom motifs or replacements, combinatorial
            concat_hairpins = TempFileHandler.write_temp_file(
                "\n".join(
                    [
                        x.csv_format
                        for x in self.motifs
                        if x.loop_type == LoopType.HAIRPIN
                    ]
                )
            )
            concat_internals = TempFileHandler.write_temp_file(
                "\n".join(
                    [
                        x.csv_format
                        for x in self.motifs
                        if x.loop_type == LoopType.INTERNAL
                    ]
                )
            )
            concat_bulges = TempFileHandler.write_temp_file(
                "\n".join(
                    [x.csv_format for x in self.motifs if x.loop_type == LoopType.BULGE]
                )
            )
            # First custom behavior, if any custom motifs or replacements are provided, combinatorial
            return [MotifCall(concat_hairpins, concat_internals, concat_bulges)]
        else:
            # Single Motif Behavior, first we sort Motifs by their abbreviations
            motif_lists: list[list[Motif]] = []
            for motif in {x.abbreviation for x in self.motifs}:
                motif_lists.append(
                    [
                        x
                        for x in self.motifs
                        if x.abbreviation == motif and x.loop_type == LoopType.HAIRPIN
                    ]
                )
                motif_lists.append(
                    [
                        x
                        for x in self.motifs
                        if x.abbreviation == motif and x.loop_type == LoopType.INTERNAL
                    ]
                )
                motif_lists.append(
                    [
                        x
                        for x in self.motifs
                        if x.abbreviation == motif and x.loop_type == LoopType.BULGE
                    ]
                )
            # Create MotifSets and filter out empty lists
            motif_sets: list[MotifSet] = [
                MotifSet.from_list(x) for x in motif_lists if len(x) > 0
            ]
            # Write Down each motif set into its own file
            for motif_set in motif_sets:
                motif_set.write_self()
            return [
                MotifHandler._split_calls(x.filepath, x.loop_type) for x in motif_sets
            ]

    @staticmethod
    def _split_calls(
        filepath: Path,
        loop_type: LoopType,
    ) -> MotifCall:
        empty_csv = (
            RNAMOTIFOLD_ROOT_DIR
            / "RNALoops"
            / "Misc"
            / "Applications"
            / "RNAmotiFold"
            / "motifs"
            / "versions"
            / "empty.csv"
        )
        match loop_type:
            case LoopType.HAIRPIN:
                return MotifCall(filepath, empty_csv, empty_csv)
            case LoopType.INTERNAL:
                return MotifCall(empty_csv, filepath, empty_csv)
            case LoopType.BULGE:
                return MotifCall(empty_csv, empty_csv, filepath)
            case LoopType.ALL:
                raise ValueError("Cant process motif of type ALL")
