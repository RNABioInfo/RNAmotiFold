import logging
from dataclasses import dataclass
from os import cpu_count
from pathlib import Path
from typing import Any, Literal, Self

from Bio import SeqRecord

import RNAmotiFold.input.cli
from RNAmotiFold.bgap_rna import (
    call_handler,
    input_handler,
    motif_handler,
    subprocess_handler,
)
from RNAmotiFold.Misc.EnumTypes import *
from RNAmotiFold.results.outputs import AlgorithmOutput
from RNAmotiFold.results.results import BaseResult


class Bunch:
    def __init__(self, d: Any = None):
        if d is not None:
            self.__dict__.update(d)


@dataclass
class Predictions:
    @classmethod
    def from_dict(cls, inputs: dict[str, Any]) -> Self:
        allowed = (
            "algorithm",
            "temperature",
            "motif_source",
            "motif_orientation",
            "output",
            "rna_3d_motif_atlas_version",
            "motif_string",
            "custom_hairpins",
            "custom_internals",
            "custom_bulges",
            "replace_hairpins",
            "replace_internals",
            "replace_bulges",
            "allowLonelyBasepairs",
            "single_motif_mode",
            "parallel_processes",
            "separator",
            "shape_level",
            "kvalue",
            "subopt",
            "subopt_energy_range_absolute",
            "subopt_energy_range_percent",
            "motif_weight",
            "motif_fraction",
            "pfc",
            "low_prob_filter",
        )
        df = {k: v for k, v in inputs.items() if k in allowed and v is not None}
        return cls(**df)

    # Necessary param
    algorithm: AlgorithmType = AlgorithmType.RNAMOTIFOLD
    ###GENERAL PARAMS####
    temperature: float = 37.0
    motif_source: MotifSource = MotifSource.BGSU
    motif_orientation: MotifOrientation = MotifOrientation.FORWARD
    output: Path | None = None
    rna_3d_motif_atlas_version: str = "4_14"
    motif_string: str = ""
    custom_hairpins: Path | None = None
    custom_internals: Path | None = None
    custom_bulges: Path | None = None
    replace_hairpins: bool = False
    replace_internals: bool = False
    replace_bulges: bool = False
    allowLonelyBasepairs: AllowLonelyBP = AllowLonelyBP.NO_LP
    single_motif_mode: bool = False
    parallel_processes: int = 1
    separator: str = ","

    ##RNAMOSH PARAMS###
    shape_level: ShapeLevel = ShapeLevel.LEVEL5

    ###MFE PARAMS###
    # Classified only
    kvalue: int = 5

    # Subopt only
    subopt: bool = False
    subopt_energy_range_absolute: float = 0.0
    subopt_energy_range_percent: float = 0.1

    # RNAmotiAlign Only
    motif_weight: float = 1.0
    motif_fraction: float = 0.7

    ###PFC PARAMS###
    pfc: bool = False
    low_prob_filter: float = 0.000001

    @property
    def prediction_type(self):
        if self.algorithm == AlgorithmType.RNAMOTIALIGN:
            return PredictionType.ALI
        elif self.pfc:
            return PredictionType.PFC
        else:
            return PredictionType.MFE


def rnamotifold(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    kvalue: int | None = 5,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTIFOLD,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        output=output,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
       allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        kvalue=kvalue,
        subopt=False,
        subopt_energy_range_absolute=None,
        subopt_energy_range_percent=None,
        motif_weight=None,
        motif_fraction=None,
        pfc=False,
        low_prob_filter=None,
    )


def rnamotifold_subopt(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    subopt_energy_range_absolute: float | None = None,
    subopt_energy_range_percent: float | None = None,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    if (
        subopt_energy_range_absolute is not None
        and subopt_energy_range_percent is not None
    ):
        raise ValueError(
            "Please only set either an absolute or a percent energy range!"
        )
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTIFOLD,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
       allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        kvalue=None,
        subopt=True,
        subopt_energy_range_absolute=subopt_energy_range_absolute,
        subopt_energy_range_percent=subopt_energy_range_percent,
        motif_weight=None,
        motif_fraction=None,
        pfc=False,
        low_prob_filter=None,
    )


def rnamotifold_pfc(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    pfc_probability_filter: float | None = None,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTIFOLD,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        kvalue=None,
        subopt=False,
        subopt_energy_range_absolute=None,
        subopt_energy_range_percent=None,
        motif_weight=None,
        motif_fraction=None,
        pfc=True,
        low_prob_filter=pfc_probability_filter,
    )


def rnamosh(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    kvalue: int | None = 5,
    shape_level:Literal[1,2,3,4,5] = 5,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOSH,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        output=output,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        kvalue=kvalue,
        subopt=False,
        subopt_energy_range_absolute=None,
        subopt_energy_range_percent=None,
        motif_weight=None,
        motif_fraction=None,
        pfc=False,
        low_prob_filter=None,
    )


def rnamosh_subopt(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    subopt_energy_range_absolute: float | None = None,
    subopt_energy_range_percent: float | None = None,
    shape_level:Literal[1,2,3,4,5] = 5,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    if (
        subopt_energy_range_absolute is not None
        and subopt_energy_range_percent is not None
    ):
        raise ValueError(
            "Please only set either an absolute or a percent energy range!"
        )
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOSH,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel(shape_level),
        kvalue=None,
        subopt=True,
        subopt_energy_range_absolute=subopt_energy_range_absolute,
        subopt_energy_range_percent=subopt_energy_range_percent,
        motif_weight=None,
        motif_fraction=None,
        pfc=False,
        low_prob_filter=None,
    )


def rnamosh_pfc(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    pfc_probability_filter: float | None = None,
    shape_level:Literal[1,2,3,4,5] = 5,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOSH,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel(shape_level),
        kvalue=None,
        subopt=False,
        subopt_energy_range_absolute=None,
        subopt_energy_range_percent=None,
        motif_weight=None,
        motif_fraction=None,
        pfc=True,
        low_prob_filter=pfc_probability_filter,
    )


def rnamotices(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    subopt_energy_range_absolute: float | None = None,
    subopt_energy_range_percent: float | None = None,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    if (
        subopt_energy_range_absolute is not None
        and subopt_energy_range_percent is not None
    ):
        raise ValueError(
            "Please only set either an absolute or a percent energy range!"
        )
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTICES,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        kvalue=None,
        subopt=True,
        subopt_energy_range_absolute=subopt_energy_range_absolute,
        subopt_energy_range_percent=subopt_energy_range_percent,
        motif_weight=None,
        motif_fraction=None,
        pfc=False,
        low_prob_filter=None,
    )


def rnamotices_subopt(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    subopt_energy_range_absolute: float | None = None,
    subopt_energy_range_percent: float | None = None,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    if (
        subopt_energy_range_absolute is not None
        and subopt_energy_range_percent is not None
    ):
        raise ValueError(
            "Please only set either an absolute or a percent energy range!"
        )
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTICES,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        kvalue=None,
        subopt=True,
        subopt_energy_range_absolute=subopt_energy_range_absolute,
        subopt_energy_range_percent=subopt_energy_range_percent,
        motif_weight=None,
        motif_fraction=None,
        pfc=False,
        low_prob_filter=None,
    )


def rnamotices_pfc(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    pfc_probability_filter: float | None = None,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    output: Path | None = None,
    separator: str | None = ",",
) -> AlgorithmOutput:
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTICES,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        kvalue=None,
        subopt=False,
        subopt_energy_range_absolute=None,
        subopt_energy_range_percent=None,
        motif_weight=None,
        motif_fraction=None,
        pfc=True,
        low_prob_filter=pfc_probability_filter,
    )


def rnamotialign(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    single_motif_mode: bool = False,
    allowLonelyBasepairs: Literal[0, 1, 2] | None = 0,
    motif_string: str = "",
    kvalue: int | None = 5,
    motif_weight: float | None = 1.0,
    motif_fraction: float | None = 0.7,
    *,
    id: str = "N/A",
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 1,
    output: Path | None = None,
    rna_3d_motif_atlas_version: str = "4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes: int | None = cpu_count(),
    separator: str | None = ",",
) -> AlgorithmOutput:
    return _predict(
        input=input,
        id=id,
        kvalue=kvalue,
        motif_fraction=motif_fraction,
        motif_weight=motif_weight,
        algorithm=AlgorithmType.RNAMOTIALIGN,
        temperature=temperature,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allowLonelyBasepairs=AllowLonelyBP(allowLonelyBasepairs),
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=ShapeLevel.LEVEL5,
        subopt=False,
        subopt_energy_range_absolute=None,
        subopt_energy_range_percent=None,
        pfc=False,
    )


def _predict(
    input: (str | SeqRecord.SeqRecord | Path | list[SeqRecord.SeqRecord] | list[str]),
    id: str = "N/A",
    algorithm: AlgorithmType = AlgorithmType.RNAMOTIFOLD,
    temperature: float = 37.0,
    motif_source: MotifSource = MotifSource.BGSU,
    motif_orientation: MotifOrientation = MotifOrientation.FORWARD,
    output: Path | None = None,
    rna_3d_motif_atlas_version: str = "4_14",
    motif_string: str = "",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    allowLonelyBasepairs: AllowLonelyBP = AllowLonelyBP.NO_LP,
    single_motif_mode: bool = False,
    parallel_processes: int | None = cpu_count(),
    separator: str | None = ",",
    ##RNAMOSH PARAMS###
    shape_level: ShapeLevel = ShapeLevel.LEVEL5,
    ###MFE PARAMS###
    # Classified only
    kvalue: int | None = 5,
    # Subopt only
    subopt: bool | None = False,
    subopt_energy_range_absolute: float | None = None,
    subopt_energy_range_percent: float | None = 5.0,
    # RNAmotiAlign Only
    motif_weight: float | None = 1.0,
    motif_fraction: float | None = 0.7,
    ###PFC PARAMS###
    pfc: bool = False,
    low_prob_filter: float | None = None,
) -> AlgorithmOutput:
    logging.basicConfig(handlers=[logging.NullHandler()])
    if parallel_processes is None:
        parallel_processes = 1
        print(
            f"Could not automatically read the number of CPU cores available, playing it save and setting to {parallel_processes}. If you want to utilize parallelized predictions, please set the number of parallel processes"
        )
    bruh = Bunch(locals())
    Predict = Predictions.from_dict(vars(bruh))
    not_installed = RNAmotiFold.input.cli._check_all_algorithms()
    if len(not_installed) > 0:
        raise FileNotFoundError(
            "Algorithms are not set up, re-run RNAmotiFold_installer"
        )
    run_process_type = GetProcessType(Predict.algorithm)
    BaseResult.separator = Predict.separator
    input_maker = input_handler.InputHandler(run_process_type, None)
    call_maker = call_handler.CallHandler(
        Predict.algorithm,
        Predict.motif_source,
        Predict.motif_orientation,
        Predict.kvalue,
        Predict.shape_level,
        Predict.subopt_energy_range_absolute,
        Predict.temperature,
        Predict.subopt_energy_range_percent,
        Predict.allowLonelyBasepairs,
        Predict.subopt,
        Predict.pfc,
        Predict.low_prob_filter,
        Predict.motif_weight,
        Predict.motif_fraction,
        Predict.rna_3d_motif_atlas_version,
    )
    main_call = call_maker.make_main_call()
    motif_call_maker = motif_handler.MotifHandler(
        Predict.motif_string,
        Predict.single_motif_mode,
        Predict.custom_hairpins,
        Predict.custom_internals,
        Predict.custom_bulges,
        Predict.replace_hairpins,
        Predict.replace_internals,
        Predict.replace_bulges,
        Predict.rna_3d_motif_atlas_version,
    )
    motif_calls = motif_call_maker.create_motif_calls()

    subprocess_manager = subprocess_handler.make_handler(
        Predict.parallel_processes,
        Predict.output,
        Predict.prediction_type,
    )
    full_calls = RNAmotiFold.input.cli.combine_calls(main_call, motif_calls)
    inputs = input_maker.script_read_input(run_process_type, input, id)
    commands = RNAmotiFold.input.cli.create_inputs(full_calls, inputs)
    results = subprocess_manager.auto_run(commands, len(motif_calls), no_print=True)
    return results
