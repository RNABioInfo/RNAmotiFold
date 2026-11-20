from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple, Self

from RNAmotiFold.bgap_rna.tempfile_handler import TempFileHandler
from RNAmotiFold.Misc.EnumTypes import LoopType, MotifOrigin
from RNAmotiFold.Misc.Helpers import flatten

from RNAmotiFold import EMPTY_CSV_PATH, RNALOOPS_PATH

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
    def from_list(cls, motlist: list[Motif]) -> Self:
        types = {x.loop_type for x in motlist}
        if len(types) > 1:
            raise ValueError("Motifs of different loop types grouped together")
        abbreviations = {x.abbreviation for x in motlist}
        if len(abbreviations) > 1:
            raise ValueError("Motifs of different abberviations grouped together")
        return cls(
            motif_list=motlist,
            abbreviation=abbreviations.pop(),
            loop_type=types.pop(),
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
class PotentialMotif:
    sequence: str
    abbreviation: str

    @classmethod
    def from_motif_file_line(cls, line: str) -> Self:
        split = line.split(",")
        if len(split) != 2:
            raise ValueError(
                f"Motif {line} is not in the correct format, should be 'sequence,abbreviation'"
            )
        return cls(split[0].strip(), split[1].strip())


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
        cls, motif_file_line: str, loop_type: LoopType, origin: MotifOrigin
    ) -> Self:
        """Creates a Motif object from a string, the string should be in the format of 'sequence,abbreviation'"""
        mot = PotentialMotif.from_motif_file_line(motif_file_line)
        match loop_type:
            case LoopType.HAIRPIN | LoopType.BULGE:
                if "$" in mot.sequence:
                    raise ValueError(
                        f"Internal sequence detected in hairpin or bulge loop file: {motif_file_line}"
                    )
            case LoopType.INTERNAL:
                if "$" not in mot.sequence:
                    raise ValueError(
                        f"Internal sequence without $ splitter detected: {motif_file_line}"
                    )
        return cls(mot.sequence, mot.abbreviation, loop_type, origin)

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
    def call(self) -> list[str]:
        returnlist: list[str] = []
        if self.hairpin_path is not None:
            returnlist.extend(["-X", f"{self.hairpin_path!s}", "-L", "1"])
        if self.internal_path is not None:
            returnlist.extend(["-Y", f"{self.internal_path!s}", "-E", "1"])
        if self.bulge_path is not None:
            returnlist.extend(["-Z", f"{self.bulge_path!s}", "-G", "1"])
        return returnlist


class ReplaceStatus(NamedTuple):
    replace_hairpins: bool
    replace_internals: bool
    replace_bulges: bool

    def any_replaced(self) -> bool:
        return any(
            x
            for x in [
                self.replace_hairpins,
                self.replace_internals,
                self.replace_bulges,
            ]
        )


class CustomMotifs(NamedTuple):
    hairpins: list[Motif]
    internals: list[Motif]
    bulges: list[Motif]


def get_custom_motifs(
    custom_path: Path | None, motif_type: LoopType
) -> list[Motif]:
    """Returns a list of Motif objects from a given file path"""
    if custom_path is None:
        return []
    return list(
        Motif.from_file(custom_path, motif_type, MotifOrigin.CUSTOM)
    )


def _get_customs(
    custom_hairpin_filepath: Path | None,
    custom_internal_path: Path | None,
    custom_bulge_path: Path | None,
) -> CustomMotifs:
    return CustomMotifs(
        get_custom_motifs(custom_hairpin_filepath, LoopType.HAIRPIN),
        get_custom_motifs(custom_internal_path, LoopType.INTERNAL),
        get_custom_motifs(custom_bulge_path, LoopType.BULGE),
    )


def _filter_default_motifs(
    potential_motifs: list[Motif], replace_bool: bool, loop_type: LoopType
):
    if not replace_bool:
        return [x for x in potential_motifs if x.loop_type == loop_type]
    else:
        return []


def create_motif_subcalls(
    hairpin_path: Path | None,
    internal_path: Path | None,
    bulge_path: Path | None,
    replace_hairpins: bool,
    replace_internals: bool,
    replace_bulges: bool,
    version: str,
    motif_string: str,
    single_motif_mode: bool,
):
    replacements = ReplaceStatus(
        replace_hairpins, replace_internals, replace_bulges
    )
    if (
        all(x is None for x in [hairpin_path, internal_path, bulge_path])
        and not replacements.any_replaced()
        and len(motif_string) == 0
        and not single_motif_mode
    ):
        return []  # instant return of no calls if no paths are set , not single mode ,no replacements and no motif string
    customs = _get_customs(hairpin_path, internal_path, bulge_path)
    default_motifs = flatten([
        list(
            Motif.from_file(
                x,
                LoopType(str(x.parent.name[:-1])),
                origin=MotifOrigin.DEFAULT,
            )
        )
        for x in get_motif_file_paths(version)
    ])
    motifs_list = customs.hairpins + customs.internals + customs.bulges
    motifs_list.extend(
        _filter_default_motifs(
            default_motifs, replacements.replace_hairpins, LoopType.HAIRPIN
        )
    )
    motifs_list.extend(
        _filter_default_motifs(
            default_motifs,
            replacements.replace_internals,
            LoopType.INTERNAL,
        )
    )
    motifs_list.extend(
        _filter_default_motifs(
            default_motifs, replacements.replace_bulges, LoopType.BULGE
        )
    )
    if len(motif_string) == 0:
        motifs = motifs_list
    else:
        motifs = [x for x in motifs_list if x.abbreviation in motif_string]
    if (
        all(x.origin == MotifOrigin.DEFAULT for x in motifs)
        and not replacements.any_replaced()
        and not single_motif_mode
        and len(motif_string) == 0
    ):
        return []  # Default behavior if motifs is either empty or they are all defaults and nothing was replaced -> No specific motif subcalls
    if single_motif_mode:
        return _single_motif_calls(motifs)
    else:
        return _combinatorial_motif_call(motifs)


def _combinatorial_motif_call(
    motifs: list[Motif],
) -> list[MotifCall]:
    # New implementation pre-filters motifs, so self.motifs only contains the motifs we want.
    # First case: Combinatorial where we can just put everything with the same motif type together

    concat_hairpins = TempFileHandler.write_temp_file(
        "\n".join([
            x.csv_format for x in motifs if x.loop_type == LoopType.HAIRPIN
        ])
    )
    concat_internals = TempFileHandler.write_temp_file(
        "\n".join([
            x.csv_format for x in motifs if x.loop_type == LoopType.INTERNAL
        ])
    )
    concat_bulges = TempFileHandler.write_temp_file(
        "\n".join([
            x.csv_format for x in motifs if x.loop_type == LoopType.BULGE
        ])
    )
    return [MotifCall(concat_hairpins, concat_internals, concat_bulges)]


def _single_motif_calls(motifs: list[Motif]) -> list[MotifCall]:
    # Single Motif Behavior, first we sort Motifs by their abbreviations
    motif_sets: list[MotifSet] = []
    for motif in {x.abbreviation for x in motifs}:
        instances_hairpins = [
            x
            for x in motifs
            if x.abbreviation == motif and x.loop_type == LoopType.HAIRPIN
        ]
        instances_internals = [
            x
            for x in motifs
            if x.abbreviation == motif and x.loop_type == LoopType.INTERNAL
        ]
        instances_bulges = [
            x
            for x in motifs
            if x.abbreviation == motif and x.loop_type == LoopType.BULGE
        ]
        motif_sets.extend([
            MotifSet.from_list(x)
            for x in [
                instances_hairpins,
                instances_internals,
                instances_bulges,
            ]
            if len(x) > 0
        ])
    # Write Down each motif set into its own file
    for motif_set in motif_sets:
        motif_set.write_self()
    return [_split_calls(x.filepath, x.loop_type) for x in motif_sets]


def _split_calls(
    filepath: Path,
    loop_type: LoopType,
) -> MotifCall:
    match loop_type:
        case LoopType.HAIRPIN:
            return MotifCall(filepath, EMPTY_CSV_PATH, EMPTY_CSV_PATH)
        case LoopType.INTERNAL:
            return MotifCall(EMPTY_CSV_PATH, filepath, EMPTY_CSV_PATH)
        case LoopType.BULGE:
            return MotifCall(EMPTY_CSV_PATH, EMPTY_CSV_PATH, filepath)
        case LoopType.ALL:
            raise ValueError("Cant process motif of type ALL")
