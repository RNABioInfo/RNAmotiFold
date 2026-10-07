from typing import Literal
from pathlib import Path
from os import cpu_count
from Bio import SeqRecord
import RNAmotiFold
import RNAmotiFold.bgap_rna
import RNAmotiFold.bgap_rna.call_handler
import RNAmotiFold.bgap_rna.input_handler
import RNAmotiFold.bgap_rna.motif_handler
import RNAmotiFold.bgap_rna.subprocess_handler
import RNAmotiFold.input
import RNAmotiFold.input.cli
import RNAmotiFold.results
import RNAmotiFold.results.base_result
from dataclasses import dataclass
from typing import Any
import logging

class Bunch(object):
    def __init__(self, d=None):
        if d is not None: self.__dict__.update(d)

@dataclass
class Predictions:
    @classmethod
    def from_dict(cls,inputs:dict[str,Any]):
        allowed = ("algorithm","temperature","motif_source","motif_orientation","output","rna_3d_motif_atlas_version","motif_string","custom_hairpins",
                   "custom_internals","custom_bulges","replace_hairpins","repalce_internals","replace_bulges","allowLonelyBasepairs","single_motif_mode",
                   "parallel_processes","separator","shape_level","kvalue","subopt","subopt_energy_range_absolute","subopt_energy_range_percent","motif_weight",
                   "motif_fraction","pfc","low_prob_filter")
        df = {k:v for k,v in inputs.items() if k in allowed and v is not None}
        return cls(**df)

    #Necessary param
    algorithm:Literal["RNAmotiFold","RNAmotiCes","RNAmotiAlign","RNAmoSh"] = "RNAmotiFold"
    ###GENERAL PARAMS####
    temperature: float = 37.0
    motif_source:Literal[1,2,3] = 1
    motif_orientation:Literal[1,2,3] = 3
    output:Path|None=None
    rna_3d_motif_atlas_version:str = "4_14"
    motif_string:str = ""
    custom_hairpins: Path | None = None
    custom_internals: Path | None = None
    custom_bulges: Path | None = None
    replace_hairpins: bool = False
    replace_internals: bool = False
    replace_bulges: bool = False
    allowLonelyBasepairs: Literal[0, 1, 2] = 0
    single_motif_mode:bool = False
    parallel_processes:int = 1
    separator:str = ","

    ##RNAMOSH PARAMS###
    shape_level:Literal[1,2,3,4,5] = 5

    ###MFE PARAMS###
    #Classified only
    kvalue:int = 5

    #Subopt only 
    subopt:bool = False
    subopt_energy_range_absolute:float = 0.0
    subopt_energy_range_percent:float = 10.0

    #RNAmotiAlign Only
    motif_weight:float = 1.0
    motif_fraction:float = 0.7

    ###PFC PARAMS###
    pfc:bool = False
    low_prob_filter: float = 0.000001


    @property
    def process_type(self):
        if self.algorithm == "RNAmotiAlign":
            return "ali"
        else:
            return "single"
        
    @property
    def prediction_type(self):
        if self.algorithm == "RNAmotiAlign":
            return "ali"
        elif self.pfc:
            return "pfc"
        else:
            return "mfe"

def rnamotifold(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    kvalue:int|None = 5,
    *,
    ID:str="N/A",
    merge_results:bool=False,
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
    return _predict(input=input,
                    ID=ID,
                    merge_results=merge_results,
                    algorithm="RNAmotiFold",
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
                    allowLonelyBasepairs=allowLonelyBasepairs,
                    single_motif_mode=single_motif_mode,
                    parallel_processes=parallel_processes,
                    separator=separator,
                    shape_level=None,
                    kvalue=kvalue,
                    subopt=False,
                    subopt_energy_range_absolute=None,
                    subopt_energy_range_percent=None,
                    motif_weight=None,
                    motif_fraction=None,
                    pfc=False,
                    low_prob_filter=None)

def rnamotifold_subopt(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    subopt_energy_range_absolute:float|None = None,
    subopt_energy_range_percent:float|None = None,
    *,
    ID:str="N/A",
    merge_results:bool=False,
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
        if subopt_energy_range_absolute is not None and subopt_energy_range_percent is not None:
            raise ValueError("Please only set either an absolute or a percent energy range!") 
        return _predict(input=input,
                        ID=ID,
                        merge_results=merge_results,
                        algorithm="RNAmotiFold",
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
                        allowLonelyBasepairs=allowLonelyBasepairs,
                        single_motif_mode=single_motif_mode,
                        parallel_processes=parallel_processes,
                        output=output,
                        separator=separator,
                        shape_level=None,
                        kvalue=None,
                        subopt=True,
                        subopt_energy_range_absolute=subopt_energy_range_absolute,
                        subopt_energy_range_percent=subopt_energy_range_percent,
                        motif_weight=None,
                        motif_fraction=None,
                        pfc=False,
                        low_prob_filter=None)

def rnamotifold_pfc(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    pfc_probability_filter: float|None = None,
    *,
    ID:str="N/A",
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
    return _predict(input=input,
                            ID=ID,
                            merge_results=False,
                            algorithm="RNAmotiFold",
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
                            allowLonelyBasepairs=allowLonelyBasepairs,
                            single_motif_mode=single_motif_mode,
                            parallel_processes=parallel_processes,
                            output=output,
                            separator=separator,
                            shape_level=None,
                            kvalue=None,
                            subopt=False,
                            subopt_energy_range_absolute=None,
                            subopt_energy_range_percent=None,
                            motif_weight=None,
                            motif_fraction=None,
                            pfc=True,
                            low_prob_filter=pfc_probability_filter)
def rnamosh(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    kvalue:int|None = 5,
    *,
    ID:str="N/A",
    merge_results:bool=False,
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
    return _predict(input=input,
                    ID=ID,
                    merge_results=merge_results,
                    algorithm="RNAmoSh",
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
                    allowLonelyBasepairs=allowLonelyBasepairs,
                    single_motif_mode=single_motif_mode,
                    parallel_processes=parallel_processes,
                    separator=separator,
                    shape_level=None,
                    kvalue=kvalue,
                    subopt=False,
                    subopt_energy_range_absolute=None,
                    subopt_energy_range_percent=None,
                    motif_weight=None,
                    motif_fraction=None,
                    pfc=False,
                    low_prob_filter=None)

def rnamosh_subopt(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    subopt_energy_range_absolute:float|None = None,
    subopt_energy_range_percent:float|None = None,
    *,
    ID:str="N/A",
    merge_results:bool=False,
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
        if subopt_energy_range_absolute is not None and subopt_energy_range_percent is not None:
            raise ValueError("Please only set either an absolute or a percent energy range!") 
        return _predict(input=input,
                        ID=ID,
                        merge_results=merge_results,
                        algorithm="RNAmoSh",
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
                        allowLonelyBasepairs=allowLonelyBasepairs,
                        single_motif_mode=single_motif_mode,
                        parallel_processes=parallel_processes,
                        output=output,
                        separator=separator,
                        shape_level=None,
                        kvalue=None,
                        subopt=True,
                        subopt_energy_range_absolute=subopt_energy_range_absolute,
                        subopt_energy_range_percent=subopt_energy_range_percent,
                        motif_weight=None,
                        motif_fraction=None,
                        pfc=False,
                        low_prob_filter=None)

def rnamosh_pfc(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    pfc_probability_filter: float|None = None,
    *,
    ID:str="N/A",
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
    return _predict(input=input,
                            ID=ID,
                            merge_results=False,
                            algorithm="RNAmoSh",
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
                            allowLonelyBasepairs=allowLonelyBasepairs,
                            single_motif_mode=single_motif_mode,
                            parallel_processes=parallel_processes,
                            output=output,
                            separator=separator,
                            shape_level=None,
                            kvalue=None,
                            subopt=False,
                            subopt_energy_range_absolute=None,
                            subopt_energy_range_percent=None,
                            motif_weight=None,
                            motif_fraction=None,
                            pfc=True,
                            low_prob_filter=pfc_probability_filter)

def rnamotices(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    subopt_energy_range_absolute:float|None = None,
    subopt_energy_range_percent:float|None = None,
    *,
    ID:str="N/A",
    merge_results:bool=False,
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
        if subopt_energy_range_absolute is not None and subopt_energy_range_percent is not None:
            raise ValueError("Please only set either an absolute or a percent energy range!") 
        return _predict(input=input,
                        ID=ID,
                        merge_results=merge_results,
                        algorithm="RNAmotiCes",
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
                        allowLonelyBasepairs=allowLonelyBasepairs,
                        single_motif_mode=single_motif_mode,
                        parallel_processes=parallel_processes,
                        output=output,
                        separator=separator,
                        shape_level=None,
                        kvalue=None,
                        subopt=True,
                        subopt_energy_range_absolute=subopt_energy_range_absolute,
                        subopt_energy_range_percent=subopt_energy_range_percent,
                        motif_weight=None,
                        motif_fraction=None,
                        pfc=False,
                        low_prob_filter=None)

def rnamotices_subopt(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    subopt_energy_range_absolute:float|None = None,
    subopt_energy_range_percent:float|None = None,
    *,
    ID:str="N/A",
    merge_results:bool=False,
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
        if subopt_energy_range_absolute is not None and subopt_energy_range_percent is not None:
            raise ValueError("Please only set either an absolute or a percent energy range!") 
        return _predict(input=input,
                        ID=ID,
                        merge_results=merge_results,
                        algorithm="RNAmotiCes",
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
                        allowLonelyBasepairs=allowLonelyBasepairs,
                        single_motif_mode=single_motif_mode,
                        parallel_processes=parallel_processes,
                        output=output,
                        separator=separator,
                        shape_level=None,
                        kvalue=None,
                        subopt=True,
                        subopt_energy_range_absolute=subopt_energy_range_absolute,
                        subopt_energy_range_percent=subopt_energy_range_percent,
                        motif_weight=None,
                        motif_fraction=None,
                        pfc=False,
                        low_prob_filter=None)

def rnamotices_pfc(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    pfc_probability_filter: float|None = None,
    *,
    ID:str="N/A",
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    output:Path|None=None,
    separator:str|None = ",",
    ):
    return _predict(input=input,
                            ID=ID,
                            merge_results=False,
                            algorithm="RNAmotiCes",
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
                            allowLonelyBasepairs=allowLonelyBasepairs,
                            single_motif_mode=single_motif_mode,
                            parallel_processes=parallel_processes,
                            output=output,
                            separator=separator,
                            shape_level=None,
                            kvalue=None,
                            subopt=False,
                            subopt_energy_range_absolute=None,
                            subopt_energy_range_percent=None,
                            motif_weight=None,
                            motif_fraction=None,
                            pfc=True,
                            low_prob_filter=pfc_probability_filter)


def rnamotialign(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    single_motif_mode:bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    motif_string:str="",
    kvalue:int|None = 5,
    motif_weight:float|None = 1.0,
    motif_fraction:float|None = 0.7,
    *,
    ID:str="N/A",
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 1,
    output:Path|None=None,
    rna_3d_motif_atlas_version:str="4_14",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    parallel_processes:int|None = cpu_count(),
    separator:str|None = ",",
    ):pass

def _predict(
    input:str|SeqRecord.SeqRecord|Path|list[SeqRecord.SeqRecord]|list[str],
    ID:str="N/A",
    merge_results:bool=False, 
    algorithm:Literal["RNAmotiFold","RNAmotiCes","RNAmotiAlign","RNAmoSh"]|None=None,
    temperature: float = 37.0,
    motif_source:Literal[1,2,3] = 1,
    motif_orientation:Literal[1,2,3] = 3,
    output:Path|None=None,
    rna_3d_motif_atlas_version:str="4_14",
    motif_string:str="",
    custom_hairpins: Path | None = None,
    custom_internals: Path | None = None,
    custom_bulges: Path | None = None,
    replace_hairpins: bool = False,
    replace_internals: bool = False,
    replace_bulges: bool = False,
    allowLonelyBasepairs:Literal[0,1,2] | None = 0,
    single_motif_mode:bool = False,
    parallel_processes:int|None = cpu_count(),
    separator:str|None = ",",

    ##RNAMOSH PARAMS###
    shape_level:Literal[1,2,3,4,5]|None = 5,

    ###MFE PARAMS###
    #Classified only
    kvalue:int|None = 5,

    #Subopt only 
    subopt:bool|None = False,
    subopt_energy_range_absolute:float|None = None,
    subopt_energy_range_percent:float|None = 5.0,

    #RNAmotiAlign Only
    motif_weight:float|None = 1.0,
    motif_fraction:float|None = 0.7,

    ###PFC PARAMS###
    pfc:bool = False,
    low_prob_filter: float|None = None
    ):
    logging.basicConfig(handlers=[logging.NullHandler()])
    if parallel_processes is None:
        parallel_processes = 1
        print(f"Could not automatically read the number of CPU cores available, playing it save and setting to {parallel_processes}. If you want to utilize parallelized predictions, please set the number of parallel processes")   
    bruh = Bunch(locals())
    Predict = Predictions.from_dict(vars(bruh))
    installed = RNAmotiFold.input.cli._check_install(Predict.algorithm)
    if not installed:
        raise FileNotFoundError("Could not find binary for chosen algorithm please run RNAmotiFold from your commandline once to install algorithm binaries")
    RNAmotiFold.results.base_result._Result._separator = Predict.separator
    input_maker = RNAmotiFold.bgap_rna.input_handler._InputHandler(Predict.process_type,None)
    call_maker = RNAmotiFold.bgap_rna.call_handler._CallHandler(Predict.algorithm,Predict.motif_source,Predict.motif_orientation,Predict.kvalue,Predict.shape_level,Predict.subopt_energy_range_absolute,Predict.temperature,Predict.subopt_energy_range_percent,Predict.allowLonelyBasepairs,Predict.subopt,Predict.pfc,Predict.low_prob_filter,Predict.single_motif_mode,Predict.motif_weight,Predict.motif_fraction,Predict.rna_3d_motif_atlas_version)
    motif_call_maker = RNAmotiFold.bgap_rna.motif_handler._MotifHandler(Predict.motif_string,Predict.single_motif_mode,Predict.custom_hairpins,Predict.custom_internals,Predict.custom_bulges,Predict.replace_hairpins,Predict.replace_internals,Predict.replace_bulges,Predict.rna_3d_motif_atlas_version)
    subprocess_manager = RNAmotiFold.bgap_rna.subprocess_handler._SubprocessHandler(Predict.parallel_processes,Predict.output,Predict.prediction_type)
    full_calls = RNAmotiFold.input.cli._combine_calls(call_maker.call,motif_call_maker.motif_calls)
    inputs = input_maker.script_read_input(Predict.process_type,input,ID)#type:ignore Raises a warning because process type is a string, but the init makes sure it can only be single or ali
    commands = RNAmotiFold.input.cli._create_inputs(full_calls,inputs)
    if len(commands ) == 1:
        return [subprocess_manager.single_run(commands[0],Predict.prediction_type,merge_results)]
    return subprocess_manager.run(commands,merge_results,motif_call_maker.call_number,no_print=True)

