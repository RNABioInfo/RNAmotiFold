from pathlib import Path




class motif_handler:
    """Motif Handler Class, this is responsible for creating the motif part of the algorithm calls.
    It should in the end provide a List of every necessary combination of algorithm calls to cover the requested motifs
    Each part should be handled individually first ?
    If the motif string is empty and no motifs are given, this shouldn't be executed at all ? Just run the algorithm normally on the requested version.
    
    First: Select the motifs we want to have as part of the computation through the motif string. If the string is empty we take all motifs.
    
    

    """

    def __init__(
        self,
        version: str,
        motifs: str,
        single_motif_mode: bool,
        custom_hairpin_filepath: Path | None,
        custom_internal_filepath: Path | None,
        custom_bulge_filepath: Path | None,
        replace_hairpins: bool,
        replace_internals: bool,
        replace_bulges: bool,
    ):
        self.version = version
        self.single_motif_mode = single_motif_mode
        self.motif_string = motifs
        self.custom_hairpins = custom_hairpin_filepath
        self.replace_hairpins = replace_hairpins
        self.custom_internals = custom_internal_filepath
        self.replace_internals = replace_internals
        self.custom_bulges = custom_bulge_filepath
        self.replace_bulges = replace_bulges


    @property
    def motif_calls(self) -> list[str]:
        #First case: Motif string is empty and all three custom motif paths are not set, single motif mode is also False --> Default case, no extra stuff necessary
        if self.motif_string == "" and all([x == None for x in [self.custom_hairpins,self.custom_internals,self.custom_bulges]]) and not self.single_motif_mode:
            return []
        separate_motif_files = self.get_motif_files()
        hairpins = motif_handler._make_file_list(self.replace_hairpins,self.custom_hairpins,"hairpins",separate_motif_files,self.motif_string)
        internals = motif_handler._make_file_list(self.replace_internals,self.custom_internals,"internals",separate_motif_files,self.motif_string)
        bulges = motif_handler._make_file_list(self.replace_bulges,self.custom_bulges,"bulges",separate_motif_files,self.motif_string)
        if self.single_motif_mode:
            
 
 
        raise NotImplementedError("Haven't implemented this part yet whoopsy")
        


    def sort_sequences(self,file_list:list[Path]):
        """Read all the files in the list into one long list with all the sequences, then sort sequences by their abbreviations, write each set to a separate temp file"""
        sequences:list[str] = []
        for file in file_list:
            with open(file,"r") as motif_file:
                sequences.extend(motif_file.readlines())
        groups = {}
        for motif_instance in sequences:
            motif = motif_instance.split(",")[-1]
            groups.setdefault(motif, []).append(motif_instance)




    def get_motif_files(self) -> list[Path]:
        motif_dir_path = (
            Path(__file__)
            .resolve()
            .parent.joinpath(
                "..",
                "..",
                "submodules",
                "RNALoops",
                "Misc",
                "Applications",
                "RNAmotiFold",
                "motifs",
                "versions",
                f"{self.version}_separated",
            )
        )
        files = Path(motif_dir_path).rglob("*.csv")
        return list(files)

    @staticmethod
    def _make_file_list(
        replace_motifs: bool,
        custom_motifs: Path|None,
        motif_type: str,
        files: list[Path],
        motif_string: str,
    ) -> list[Path]:
        """Give the function the replace bool, custom motif Path (optional) and which type of motifs its working with"""
        if not replace_motifs:  # This catches both cases, 0 for no replacing and None (which also means no replacing)
            motif_paths = [x for x in files if motif_type in x.parent.name and motif_handler.check_abb(x, motif_string)]
            if custom_motifs is not None:
                motif_paths.append(custom_motifs)
            return motif_paths
        else:
            if custom_motifs is not None:
                return [custom_motifs]
            raise ValueError(
                f"Replacement of {motif_type} motifs was set to True but no custom motif file was provided, please provide a custom motif file or set replacement to False"
            )  # I kinda want to just let this slide and make the software just use no motifs in this case

    @staticmethod
    def check_abb(filepath: Path, motif_string: str):
        if len(motif_string) == 0:
            return True
        with open(filepath, "r") as motif_seq_file:
            lines = motif_seq_file.readlines()
            abbreviations = set([x.split(",")[1].strip() for x in lines])
            return all([x in motif_string for x in abbreviations])


if __name__ == "__main__":
    newhandler1 = motif_handler("4_10","",False,None,None,None,False,False,False)
    newhandler2 = motif_handler("4_10","",False,Path("/home/ubuntu/ztest/test.csv").resolve(),None,None,False,False,False)
    newhandler3 = motif_handler("4_10","",False,Path("/home/ubuntu/ztest/test.csv").resolve(),None,None,True,False,False)
    bruh1 = newhandler1.motif_calls
    bruh2 = newhandler2.motif_calls
    bruh3 = newhandler3.motif_calls