from pathlib import Path
import tempfile
import glob
import os
from typing import Literal
from src.input.parameters import ScriptParameters


class motif_handler:
    """Motif Handler Class, this is responsible for creating the motif part of the algorithm calls.
    It should in the end provide a List of every necessary combination of algorithm calls to cover the requested motifs
    Each part should be handled individually first ?
    If the motif string is empty and no motifs are given, returns am empty list.
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
        self._tmp_folders: None|tuple[Path,Path,Path] = None

    @classmethod
    def from_script_parameters(cls,params:ScriptParameters) -> "motif_handler":
        return cls(params.version,params.motif_list,params.fast_mode,params.custom_hairpins,params.custom_internals,params.custom_bulges,params.replace_hairpins,params.replace_internals,params.replace_bulges)


    @property
    def motif_calls(self) -> list[str]:
        # First case: Motif string is empty and all three custom motif paths are not set, single motif mode is also False --> Default case, no extra stuff necessary
        if (
            self.motif_string == ""
            and all(
                [
                    x == None
                    for x in [
                        self.custom_hairpins,
                        self.custom_internals,
                        self.custom_bulges,
                    ]
                ]
            )
            and not self.single_motif_mode
        ):
            return []
        separate_motif_files = self.get_motif_files()
        hairpins = motif_handler._make_file_list(
            self.replace_hairpins,
            self.custom_hairpins,
            "hairpins",
            separate_motif_files,
            self.motif_string,
        )
        internals = motif_handler._make_file_list(
            self.replace_internals,
            self.custom_internals,
            "internals",
            separate_motif_files,
            self.motif_string,
        )
        bulges = motif_handler._make_file_list(
            self.replace_bulges,
            self.custom_bulges,
            "bulges",
            separate_motif_files,
            self.motif_string,
        )
        if self.single_motif_mode:
            temp_folders = (
                motif_handler._split_sequences(
                    self.motif_string, hairpins
                ),
                motif_handler._split_sequences(
                    self.motif_string, internals
                ),
                motif_handler._split_sequences(
                    self.motif_string, bulges
                ),
            )
            self._tmp_folders = temp_folders
            # Jetzt: alle drei Folder globben, dann hab ich die Paths zu jedem einzelnen Motif separat. Danach einfach kombinieren jedes file mit 2x empty csv und die kombinationen returnen
            hairpin_calls = motif_handler._split_calls(
                glob.glob(str(temp_folders[0] / "*.tmp")), "hairpin"
            )
            internal_calls = motif_handler._split_calls(
                glob.glob(str(temp_folders[1] / "*.tmp")), "internal"
            )
            bulge_calls = motif_handler._split_calls(
                glob.glob(str(temp_folders[2] / "*.tmp")), "bulge"
            )
            return hairpin_calls + internal_calls + bulge_calls
        else:
            concat_hairpins = self._filter_concat(
                hairpins, self.motif_string
            )
            concat_internals = self._filter_concat(
                internals, self.motif_string
            )
            concat_bulges = self._filter_concat(
                bulges, self.motif_string
            )
            self._tmp_folders = (concat_hairpins.parent,concat_internals.parent,concat_bulges.parent)
            return [
                f"-X {concat_hairpins} -Y {concat_internals} -Z {concat_bulges} -L 1 -E 1 -G 1"
            ]

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
    def _split_calls(
        file_list, loop_type: Literal["hairpin", "internal", "bulge"]
    ):
        empty_csv = (
            Path(__file__).resolve().parents[2]
            / "submodules"
            / "RNALoops"
            / "Misc"
            / "Applications"
            / "RNAmotiFold"
            / "motifs"
            / "versions"
            / "empty.csv"
        )
        match loop_type:
            case "hairpin":
                return [
                    f"-X {x} -Y {empty_csv} -Z {empty_csv} -L 1 -E 1 -G 1"
                    for x in file_list
                ]
            case "internal":
                return [
                    f"-X {empty_csv} -Y {x} -Z {empty_csv} -L 1 -E 1 -G 1"
                    for x in file_list
                ]
            case "bulge":
                return [
                    f"-X {empty_csv} -Y {empty_csv} -Z {x} -L 1 -E 1 -G 1"
                    for x in file_list
                ]

    @staticmethod
    def _filter_concat(file_list, motif_string) -> Path:
        """Takes a list of motif files and a motif string, reads all the files, filters out only those in the motif string and puts them back together. If the motif string is empty it takes all sequences
        Returns the path to the new temp file with the sequences in it.
        """
        all_sequences = motif_handler._sort_sequences(
            file_list, motif_string
        )
        contents = []
        for key in all_sequences:
            contents.extend(all_sequences[key])
        temp = tempfile.NamedTemporaryFile(
            delete_on_close=False, suffix=".tmp"
        )
        temp.write("".join(contents).encode())
        return Path(temp.name).resolve()

    @staticmethod
    def _sort_sequences(file_list, motif_string):
        """Function for both single motif mode and combinatorial, reads all the files given to it, checks if abbreviations are in the motif string, and returns dictionary of
        abbreviation to all sequence variants (still including the abbreviation and the newline!). If motif string is empty, everything is kept
        """
        sequences: list[str] = []
        for file in file_list:
            with open(file, "r") as motif_file:
                sequences.extend(motif_file.readlines())
        groups: dict[str, list[str]] = {}
        for motif_instance in sequences:
            motif = motif_instance.split(",")[-1].strip()
            groups.setdefault(motif, []).append(
                motif_instance
            )  # Set default either sets the given default or returns the value if the key already exsist
        if len(motif_string) != 0:
            deletes = []
            for key in groups.keys():
                if key not in motif_string:
                    deletes.append(key)
            for char in deletes:
                del groups[char]
        return groups

    @staticmethod
    def _split_sequences(motif_string, file_list: list[Path]) -> Path:
        """Read all the files in the list into one long list with all the sequences, then sort sequences by their abbreviations,
        write each set to a separate temp file and finally return the path to a tempdir where all the tempfiles have been written.
        Remember to clean up the tempdir after to remvoe all the temp files!
        """
        tempdir: tempfile.TemporaryDirectory[str] = (
            tempfile.TemporaryDirectory(delete=False)
        )
        groups: dict[str, list[str]] = motif_handler._sort_sequences(
            file_list, motif_string
        )
        for key in groups.keys():
            motifs = "".join(groups[key])
            motif_temp = tempfile.NamedTemporaryFile(
                dir=tempdir.name,
                delete=False,
                delete_on_close=False,
                suffix=".tmp",
            )
            motif_temp.write(motifs.encode())
        return Path(tempdir.name).resolve()

    @staticmethod
    def _make_file_list(
        replace_motifs: bool,
        custom_motifs: Path | None,
        motif_type: str,
        files: list[Path],
        motif_string: str,
    ) -> list[Path]:
        """Give the function the replace bool, custom motif Path (optional) and which type of motifs its working with"""
        if (
            not replace_motifs
        ):  # This catches both cases, 0 for no replacing and None (which also means no replacing)
            motif_paths = [
                x
                for x in files
                if motif_type in x.parent.name
                and motif_handler.check_abb(x, motif_string)
            ]
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
            abbreviations = set(
                [x.split(",")[1].strip() for x in lines]
            )
            return all([x in motif_string for x in abbreviations])

    def cleanup_temp_files(self):
        if self._tmp_folders is not None:
            for folder in self._tmp_folders:
                files = glob.glob(str(folder / "*.tmp"))
                for file in files:
                    os.remove(file)