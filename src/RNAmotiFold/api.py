from decimal import Decimal
from os import cpu_count
from pathlib import Path
from typing import Literal

from Bio import SeqRecord
from Bio.Align import MultipleSeqAlignment

import RNAmotiFold.input.cli
from RNAmotiFold.bgap_rna import subprocess_handler
from RNAmotiFold.bgap_rna.call_handler import (
    main_call_from_script_params,
)
from RNAmotiFold.bgap_rna.input_handler import get_input
from RNAmotiFold.bgap_rna.motif_handler import (
    MotifCall,
    create_motif_subcalls,
)
from RNAmotiFold.bgap_rna.tempfile_handler import TempFileHandler
from RNAmotiFold.input.parameters import ScriptParameters
from RNAmotiFold.Misc.EnumTypes import *
from RNAmotiFold.Misc.ValueObjects import (
    SubOptEnergy,
    fraction,
    limited_temperature,
    pos_decimal,
    pos_int,
)
from RNAmotiFold.results.outputs import AlgorithmOutput
from RNAmotiFold.results.results import BaseResult


def rnamotifold(
    input: str | SeqRecord.SeqRecord | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    kvalue: int = 5,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a rnamotifold prediction. Arguments:
    input -- Input for RNAmotiFold, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    kvalue -- Set the number of classes to be calculated for each prediction
    Keyword Only Arguments:
    parallel_processes -- RNAmotiFold is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """

    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTIFOLD,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        output=output,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        separator=separator,
        kvalue=kvalue,
    )


def rnamotifold_subopt(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    subopt_energy: float = 5.0,
    subopt_mode: Literal["A", "P"] = "P",
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmotiFold subopt prediction. This always classifies ALL possible motif combinations and may take very long for big inputs. Arguments:
    input -- Input for RNAmotiFold, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    subopt_energy -- Specific the energy range in which to calculate all suboptimal structures, applied to each class separately. READ subopt_mode TOO!
    subopt_mode -- Set your subopt mode, either to mode 'A' for absolute or 'P' for percentage. If 'A' is set, subopt_energy is interpreted as kcal/mol, if 'P' is set subopt_energy is interpreted as percent!
    Keyword Only Arguments:
    parallel_processes -- RNAmotiFold is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTIFOLD,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        subopt=True,
        subopt_energy=subopt_energy,
        subopt_mode=subopt_mode,
    )


def rnamotifold_pfc(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    pfc_probability_filter: float = 0.01,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmotiFold pfc prediction. This always classifies ALL possible motif combinations and may take very long for big inputs. Arguments:
    input -- Input for RNAmotiFold, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    pfc_probability_filter -- Set a probability filter for classes, this filters out classes with a probability lower than the filter.
    Keyword Only Arguments:
    parallel_processes -- RNAmotiFold is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTIFOLD,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        pfc=True,
        low_prob_filter=pfc_probability_filter,
    )


def rnamosh(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    kvalue: int = 5,
    shape_level: Literal[1, 2, 3, 4, 5] = 5,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmoSh prediction. Arguments:
    input -- Input for RNAmoSh, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    kvalue -- Set the number of classes to be calculated for each prediction
    shape_level -- Set abstraction level of shapes
    Keyword Only Arguments:
    parallel_processes -- RNAmoSh is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOSH,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        output=output,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        separator=separator,
        shape_level=shape_level,
        kvalue=kvalue,
    )


def rnamosh_subopt(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    subopt_energy: float = 5.0,
    subopt_mode: Literal["A", "P"] = "P",
    shape_level: Literal[1, 2, 3, 4, 5] = 5,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmoSh subopt prediction. This always classifies ALL possible motif combinations and may take very long for big inputs. Arguments:
    input -- Input for RNAmoSh, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    shape_level -- Set abstraction level of shapes
    subopt_energy -- Specific the energy range in which to calculate all suboptimal structures, applied to each class separately. READ subopt_mode TOO!
    subopt_mode -- Set your subopt mode, either to mode 'A' for absolute or 'P' for percentage. If 'A' is set, subopt_energy is interpreted as kcal/mol, if 'P' is set subopt_energy is interpreted as percent!
    Keyword Only Arguments:
    parallel_processes -- RNAmoSh is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOSH,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=shape_level,
        subopt=True,
        subopt_energy=subopt_energy,
        subopt_mode=subopt_mode,
    )


def rnamosh_pfc(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    pfc_probability_filter: float = 0.01,
    shape_level: Literal[1, 2, 3, 4, 5] = 5,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmoSh pfc prediction. This always classifies ALL possible motif combinations and may take very long for big inputs. Arguments:
    input -- Input for RNAmoSh, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    pfc_probability_filter -- Set a probability filter for classes, this filters out classes with a probability lower than the filter.
    shape_level -- Set abstraction level of shapes
    Keyword Only Arguments:
    parallel_processes -- RNAmoSh is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOSH,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=shape_level,
        pfc=True,
        low_prob_filter=pfc_probability_filter,
    )


def rnamotices(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    kvalue: int = 5,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmotiCes prediction. Arguments:
    input -- Input for RNAmotiCes, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    kvalue -- Set the number of classes to be calculated for each prediction
    Keyword Only Arguments:
    parallel_processes -- RNAmotiCes is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTICES,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=5,
        kvalue=kvalue,
    )


def rnamotices_subopt(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    subopt_energy: float = 5.0,
    subopt_mode: Literal["A", "P"] = "P",
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmotiCes subopt prediction. This always classifies ALL possible motif combinations and may take very long for big inputs. Arguments:
    input -- Input for RNAmotiCes, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    subopt_energy -- Specific the energy range in which to calculate all suboptimal structures, applied to each class separately. READ subopt_mode TOO!
    subopt_mode -- Set your subopt mode, either to mode 'A' for absolute or 'P' for percentage. If 'A' is set, subopt_energy is interpreted as kcal/mol, if 'P' is set subopt_energy is interpreted as percent!
    Keyword Only Arguments:
    parallel_processes -- RNAmotiCes is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTICES,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        shape_level=5,
        subopt=True,
        subopt_energy=subopt_energy,
        subopt_mode=subopt_mode,
    )


def rnamotices_pfc(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    pfc_probability_filter: float = 0.01,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Run a RNAmotiCes pfc prediction. This always classifies ALL possible motif combinations and may take very long for big inputs. Arguments:
    input -- Input for RNAmotiCes, can be a string for a single sequence, SeqRecord, or a pathlib.Path for a file or folder
    single_motif_mode -- Activates single motif mode for this input
    allow_lonely_bp -- Set lone base pair mode: 0 = No lone bp, 1 = All lone bp, 2 = Allow lone bp around motifs only
    motif_string -- Sets the motifs to be used in your prediction
    pfc_probability_filter -- Set a probability filter for classes, this filters out classes with a probability lower than the filter.
    Keyword Only Arguments:
    parallel_processes -- RNAmotiCes is multiprocessed using python, set the number of parallel processes you want to run (Max your number of available CPU cores/hyperthreads)
    id -- For Single sequences (String input), set an ID for the output
    temperature -- Scale NN energy parameters to different temperature
    motif_souce -- Set motif source, 1 == BGSU RNA 3D Motif Atlas, 2 == Rfam, 3 == Both
    motif_orientation -- Set direction of motif sequences, 1 == 5' -> 3', 2 == 3' -> 5', 3 == Both
    rna3d_motif_atlas_version -- Set the RNA3D Motif Atlas version
    custom_hairpins -- Set Path to file with custom hairpin loop motif sequences
    custom_internals -- Set Path to file with custom internal loop motif sequences
    custom_bulges -- Set Path to file with custom bulge loop motif sequences
    replace_hairpins -- Set to True to replace all hairpin sequences with your custom motifs from custom_hairpins
    replace_internals -- Set to True to replace all internal sequences with your custom motifs from custom_internals
    replace_bulges -- Set to True to replace all bulge sequences with your custom motifs from custom_bulges
    output -- Set an output path if you want to directly write your results to a specific file in CSV format
    separator -- Set a separator for the CSV formatted output
    """
    return _predict(
        input=input,
        id=id,
        algorithm=AlgorithmType.RNAMOTICES,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
        pfc=True,
        low_prob_filter=pfc_probability_filter,
    )


def rnamotialign(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    single_motif_mode: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    motif_string: str = "",
    kvalue: int = 5,
    motif_weight: float = 1.0,
    motif_fraction: float = 0.7,
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
    separator: str = ",",
) -> AlgorithmOutput:
    """Work in progress, no peeking!"""
    return _predict(
        input=input,
        id=id,
        kvalue=kvalue,
        motif_fraction=motif_fraction,
        motif_weight=motif_weight,
        algorithm=AlgorithmType.RNAMOTIALIGN,
        temperature=temperature,
        motif_source=motif_source,
        motif_orientation=motif_orientation,
        rna_3d_motif_atlas_version=rna_3d_motif_atlas_version,
        motif_string=motif_string,
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        allow_lone_bp=allow_lone_bp,
        single_motif_mode=single_motif_mode,
        parallel_processes=parallel_processes,
        output=output,
        separator=separator,
    )


def _predict(
    input: str | SeqRecord.SeqRecord | MultipleSeqAlignment | Path,
    id: str = "N/A",
    algorithm: AlgorithmType = AlgorithmType.RNAMOTIFOLD,
    temperature: float = 37.0,
    motif_source: Literal[1, 2, 3] = 1,
    motif_orientation: Literal[1, 2, 3] = 3,
    output: Path | None = None,
    rna_3d_motif_atlas_version: str = "4_14",
    motif_string: str = "",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    allow_lone_bp: Literal[0, 1, 2] = 0,
    single_motif_mode: bool = False,
    parallel_processes: int | None = cpu_count(),
    separator: str = ",",
    ##RNAMOSH PARAMS###
    shape_level: Literal[1, 2, 3, 4, 5] = 5,
    ###MFE PARAMS###
    # Classified only
    kvalue: int = 5,
    # Subopt only
    subopt: bool = False,
    subopt_energy: float = 5.0,
    subopt_mode: Literal["A", "P"] = "P",
    # RNAmotiAlign Only
    motif_weight: float = 1.0,
    motif_fraction: float = 0.7,
    ###PFC PARAMS###
    pfc: bool = False,
    low_prob_filter: float = 0.001,
) -> AlgorithmOutput:
    if parallel_processes is None:
        parallel_processes = 1
        print(
            f"Could not automatically read the number of CPU cores available, playing it save and setting to {parallel_processes}. If you want to utilize parallelized predictions, please set the number of parallel processes"
        )

    rt_args = ScriptParameters(
        id=id,
        input=input,
        output=output,
        algorithm=algorithm,
        subopt=subopt,
        motif_source=MotifSource(motif_source),
        motif_orientation=MotifOrientation(motif_orientation),
        subopt_opts=SubOptEnergy(pos_decimal(subopt_energy), SubOptMode(subopt_mode)),
        kvalue=pos_int(kvalue),
        shape_level=ShapeLevel(shape_level),
        pfc=pfc,
        low_prob_filter=fraction(low_prob_filter),
        custom_hairpins=custom_hairpins,
        custom_internals=custom_internals,
        custom_bulges=custom_bulges,
        replace_hairpins=replace_hairpins,
        replace_internals=replace_internals,
        replace_bulges=replace_bulges,
        loglevel="DEBUG",
        logfile=None,
        workers=pos_int(parallel_processes),
        separator=separator,
        version=rna_3d_motif_atlas_version,
        single_motif_mode=single_motif_mode,
        motif_list=motif_string,
        motif_weight=Decimal(motif_weight),
        motif_fraction=fraction(motif_fraction),
        temperature=limited_temperature(temperature),
        basepairs=AllowLonelyBP(allow_lone_bp),
    )

    BaseResult.separator = rt_args.separator

    # Check if everything is set up properly
    not_installed = RNAmotiFold.input.cli._check_all_algorithms()
    if len(not_installed) > 0:
        raise FileNotFoundError(
            "Algorithms are not set up, re-run RNAmotiFold_installer"
        )
    main_call = main_call_from_script_params(rt_args)
    motif_calls: list[MotifCall] = create_motif_subcalls(
        rt_args.custom_hairpins,
        rt_args.custom_internals,
        rt_args.custom_bulges,
        rt_args.replace_hairpins,
        rt_args.replace_internals,
        rt_args.replace_bulges,
        rt_args.version,
        rt_args.motif_list,
        rt_args.single_motif_mode,
    )

    inputs = get_input(main_call, motif_calls, rt_args.process_type, input, id)

    subprocess_manager = subprocess_handler.make_handler(
        rt_args.workers,
        rt_args.output,
        rt_args.alg_type,
    )

    results = subprocess_manager.auto_run(inputs, no_print=True)
    TempFileHandler.cleanup()
    return results
