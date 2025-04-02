import json
from pathlib import Path
import requests
from os import scandir, remove
from dataclasses import dataclass
from collections import namedtuple
import re
import json

non_listed_conversions = {"MAD": "A"}
annotations = namedtuple("annotations", ["motif", "alignment", "annotations"])
conversion_json_path = (
    Path(__file__).resolve().parent.joinpath("data", "nucleotide_conversion.json")
)
motifs_folder_path = Path(__file__).resolve().parent.joinpath("data", "motifs")
duplicates_json_path = Path(__file__).resolve().parent.joinpath("data", "duplicates.json")

import logging

logger = logging.getLogger("update_motifs")


@dataclass
class MotifSequence:
    sequence: str
    abbreviation: str

    def __str__(self):
        return ",".join([self.sequence, self.abbreviation])


class rna3d_motif:

    rna3datlas_api_call = "http://rna.bgsu.edu/correspondence/pairwise_interactions_single?selection_type=loop_id&selection="

    @staticmethod
    def load_motif_json():
        """load motifs.json file for curated motif collection"""
        json_path = Path.joinpath(Path(__file__).resolve().parents[0], "data", "motifs.json")
        with open(json_path) as json_file:
            motifs_json = json.load(json_file, object_hook=rna3d_motif)  # type:list[rna3d_motif]
        return motifs_json

    @staticmethod
    def api_call(call, attempts: int = 6) -> requests.Response:
        i = 0
        while i < attempts:
            response = requests.get(call)
            if response.status_code == 200:
                return response
            else:
                i += 1
        raise ConnectionError(f"Could not retrieve data in {str(i)} attempts")

    @staticmethod
    def get_nucleotide_element(nucleotide: str, number: int) -> str:
        """Nucleotide element grabber that also converts non standard nucleotides using the nucleotide_conversion.json"""
        split = nucleotide.split("|")
        element = split[number]
        if number == 3:  # If you are grabbing the base from the nucleotide string
            if element not in ["G", "C", "U", "A"]:
                element = rna3d_motif.nucleotide_conversion(element)
        return element

    @staticmethod
    def nucleotide_conversion(nucleotide):
        if nucleotide in non_listed_conversions.keys():
            return non_listed_conversions[nucleotide]
        with open(conversion_json_path, "r") as file:
            conversion_json = json.load(file)
        return conversion_json[nucleotide]["standard_base"][0]

    @staticmethod
    def reverser(input_list):
        rev_list = []
        for entry in input_list:
            rev_list.append(entry[::-1])
        return rev_list

    @staticmethod
    def get_break(nucleotide_list: list) -> int:
        numbers = [rna3d_motif.get_nucleotide_element(x, 4) for x in nucleotide_list]
        chains = [rna3d_motif.get_nucleotide_element(x, 2) for x in nucleotide_list]
        for i in range(len(chains) - 1):
            if chains[i + 1] != chains[i] or abs(int(numbers[i + 1]) - int(numbers[i])) > 5:
                return i
            else:
                pass
        return -1

    @staticmethod
    def FUSION(a: str, b: str) -> str:
        """Either appends a with $ and b or just returns the one that isn't len 0"""
        if len(a) > 0 and len(b) > 0:
            sequence = a + "$" + b
        else:
            sequence = a + b
        return sequence

    @property
    def instance_ids(self):
        return flatten([[*item["alignment"]] for item in self.motifs])

    @property
    def rna3d(self):
        if hasattr(self, "_rna3d"):
            return self._rna3d
        else:
            return []

    @rna3d.setter
    def rna3d(self, new_data):
        if hasattr(self, "_rna3d"):
            self._rna3d.append(new_data)
        else:
            self._rna3d = []
            self._rna3d.append(new_data)

    @property
    def bulge_sequences(self):
        if self.loop_type == "internal":
            return [
                MotifSequence(sequence=x, abbreviation=self.abbreviation)
                for x in set(self.rna3d)
                if "$" not in x
            ]
        else:
            raise LookupError(
                "This motif does not have bulge sequences because it is not an internal loop"
            )

    @property
    def internal_sequences(self):
        if self.loop_type == "internal":
            return [
                MotifSequence(sequence=x, abbreviation=self.abbreviation)
                for x in set(self.rna3d)
                if "$" in x
            ]
        else:
            raise LookupError(
                "This motif does not have internal sequences because it is not an internal loop"
            )

    @property
    def seq_abb_tpls(self):
        return [MotifSequence(sequence=x, abbreviation=self.abbreviation) for x in set(self.rna3d)]

    @staticmethod
    def load_current_rna3datlas():
        current_hl_json = json.loads(
            rna3d_motif.api_call(
                "http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json"
            ).content.decode()
        )
        current_hl_annotations = rna3d_motif._get_annotations(
            rna3d_motif.api_call(
                "http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/annotations"
            ).content.decode()
        )
        current_il_json = json.loads(
            rna3d_motif.api_call(
                "http://rna.bgsu.edu/rna3dhub/motifs/release/il/current/json"
            ).content.decode()
        )
        current_il_annotations = rna3d_motif._get_annotations(
            rna3d_motif.api_call(
                "http://rna.bgsu.edu/rna3dhub/motifs/release/il/current/annotations"
            ).content.decode()
        )
        return (current_hl_json, current_hl_annotations, current_il_json, current_il_annotations)

    @staticmethod
    def _get_annotations(annotation_csv: str) -> list[str]:
        return_list = []
        for line in annotation_csv.split("\n"):
            split = line.split("\t")
            # Get annotations
            annotated = [x for x in split[2:] if x != ""]
            # Only keep obj if there are any annotations [not annotated means len(list) is 0, pythonic way to check list fill status]
            if annotated:
                return_list.append(
                    annotations(motif=split[0], alignment=split[1], annotations=annotated)
                )
            else:
                pass
        return return_list

    @property
    def include_json(self):
        return self._include_json

    @include_json.setter
    def include_json(self, motif_json_val):
        if motif_json_val in ["True", "true", "yes", "y"]:
            self._include_json = True
        elif motif_json_val in ["False", "false", "no", "n"]:
            self._include_json = False
        else:
            raise ValueError(
                "Include bulged not set to True or False, please decide to include sequences from the downloadable json or not."
            )

    def __init__(self, motif_json: dict):
        self.motif_name = motif_json["motif_name"]
        self.abbreviation = motif_json["abbreviation"]
        self.loop_type = motif_json["loop_type"]
        self.include_json = motif_json["include_unbulged"]
        self.annotations = []  # type:list[annotations]

    def motifs2csv(self, mot_list):
        return "\n".join([x + f",{self.abbreviation}" for x in set(mot_list)]) + "\n"

    # Takes a list of RNA3DAtlas nucleotide sequences and extracts sequences
    def sequence_from_nucleotide_list(self, alignment_nucleotide_list):
        match self.loop_type:
            case "hairpin":
                for alignment in alignment_nucleotide_list:
                    if len(alignment) - 2 > 3:
                        self.rna3d = "".join(
                            rna3d_motif.get_nucleotide_element(x, 3) for x in alignment[1:-1]
                        )
                    else:
                        pass
            case "internal":
                for alignment in alignment_nucleotide_list:
                    sequence_break = rna3d_motif.get_break(alignment)
                    # Sequence break is set to -1 if no seq break was found, leading to the sequence being discarded
                    if sequence_break > 0:
                        front = alignment[1:sequence_break]
                        back = alignment[sequence_break + 2 : -1]
                        self.rna3d = rna3d_motif.FUSION(
                            "".join(rna3d_motif.get_nucleotide_element(x, 3) for x in front),
                            "".join(rna3d_motif.get_nucleotide_element(x, 3) for x in back),
                        )
                    else:
                        pass
                        # print(alignment) #Debugging print statement

    def get_rna3d_json_sequences(self, json_file):  # Extracts sequences from the json
        collection_list = []
        for entry in json_file:
            if entry["motif_id"] in [x.motif for x in self.annotations]:
                for json_alignment in entry["alignment"]:
                    if json_alignment in [x.alignment for x in self.annotations]:
                        collection_list.append(entry["alignment"][json_alignment])
        self.sequence_from_nucleotide_list(collection_list)

    def get_rna3d_api_sequences(self):
        collection_list = []
        for entry in self.annotations:
            api_data = (
                rna3d_motif.api_call(self.rna3datlas_api_call + entry.alignment)
                .content.decode()
                .split()
            )
            nucleotides = [x for x in api_data if "|" in x]
            collection_list.append(nucleotides)
        self.sequence_from_nucleotide_list(collection_list)

    @staticmethod  # Main function for creating the hexdumb style data
    def sequences2header(seq_set: list, name: str, hexdumb_file: Path) -> str:
        joined_seq_set = "\n".join(seq_set)
        out = []
        out.append("static char {var_name}[] = {{".format(var_name=name))
        data = [joined_seq_set[i : i + 12] for i in range(0, len(joined_seq_set), 12)]
        for i, x in enumerate(data):
            line = ", ".join((["0x{val:02x}".format(val=ord(c)) for c in x]))
            out.append(
                "  {lined}{comma}".format(lined=line, comma="," if i < len(data) - 1 else "")
            )
        out.append("};")
        out.append(
            "static unsigned int {var_name}_len = {data_len};\n".format(
                var_name=name, data_len=len(joined_seq_set)
            )
        )
        with open(hexdumb_file, "a+") as writefile:
            writefile.write("\n".join(out) + "\n")

    @staticmethod  # collects the rfam_hairpins and rfam_internals, these are pre-written so it just has to read the csv files. Both are turned into lists made up of (sequence, motif abbreviation) tuples.
    def get_rfam_sequences() -> tuple[list[MotifSequence], list[MotifSequence]]:
        rfam_hairpins = motifs_folder_path.joinpath("rfam_hairpins.csv")
        rfam_internals = motifs_folder_path.joinpath("rfam_internals.csv")
        with open(rfam_hairpins) as file:
            hairpins_read = file.readlines()  # type:list[str]
        hairpins = [
            MotifSequence(
                sequence=hairpin.split(",")[0], abbreviation=hairpin.strip().split(",")[1]
            )
            for hairpin in hairpins_read
        ]
        with open(rfam_internals) as file2:
            internals_read = file2.readlines()
        internals = [
            MotifSequence(
                sequence=internal.split(",")[0], abbreviation=internal.strip().split(",")[1]
            )
            for internal in internals_read
        ]
        return (hairpins, internals)

    def get_alignments(self, annotations_list: list[annotations]):
        """This function alters the annotations_list and deletes the annotation objects that it consumed"""
        indices = []
        for idx, element in enumerate(annotations_list):
            for slot in element.annotations:
                if len(slot) >= len(self.motif_name):
                    searchd = re.search(pattern=self.motif_name.lower(), string=slot.lower())
                    if searchd is not None:
                        if any(
                            [
                                re.search(pattern=x.lower(), string=slot.lower())
                                for x in ["mini", "variation", "related", "reverse"]
                            ]
                        ):
                            pass
                        else:
                            self.annotations.append(element)
                            indices.append(idx)
                            break
                else:
                    pass
        return [j for i, j in enumerate(annotations_list) if i not in indices]


# Function that takes an input list of MotifSequence Objects writes sequences to csv. Ignores same motif duplicate sequences.
def write_csv(input_list: list[MotifSequence], filename: str):
    seen = set([str(MotSeq) for MotSeq in input_list])
    with open(motifs_folder_path.joinpath(filename), mode="w+") as file:
        for MotSeq in seen:
            file.write(MotSeq + "\n")
    return True


def flatten(xss):  # Makes a list of lists into a single list with all the items
    return [x for xs in xss for x in xs]


# Takes a list of MotifSequence Objects and creates a new list of MotifSequence Objects with reversed sequences
def sequence_reverser(seq_abb_tpls: list[MotifSequence]) -> list[MotifSequence]:
    return [
        MotifSequence(sequence=x.sequence[::-1], abbreviation=x.abbreviation) for x in seq_abb_tpls
    ]


def update_hexdumbs():  # Function iterates of motifs csv files and writes all their hexdums to submodules/RNALoops/Extensions/mot_header.hh
    """main function for updating hexdumbs. Can be used on it's own to update the mot_header.hh file in submodules/RNALoops/Extensions."""
    if Path.is_file(
        Path(__file__)
        .resolve()
        .parent.joinpath("..", "submodules", "RNALoops", "Extensions", "mot_header.hh")
    ):
        remove(
            Path(__file__)
            .resolve()
            .parent.joinpath("..", "submodules", "RNALoops", "Extensions", "mot_header.hh")
        )
    onlyFiles = [Path(f) for f in scandir(motifs_folder_path)]
    for file in onlyFiles:
        if file.suffix == ".csv":
            with open(file, "r") as file2:
                rna3d_motif.sequences2header(
                    file2.readlines(),
                    file.stem,
                    Path(__file__)
                    .resolve()
                    .parent.joinpath("..", "submodules", "RNALoops", "Extensions", "mot_header.hh"),
                )


def main():
    """Updates motif sequence collection and the hexdumb file RNAmotiFold/submodules/RNALoops/Extensions/mot_header.hh"""
    # rfam_sequences_tuples[0] = rfam_hairpins, rfam_sequences_tuples[1] = rfam_internals
    rfam_sequences_tuples = rna3d_motif.get_rfam_sequences()
    motifs = rna3d_motif.load_motif_json()
    hl_json, hl_annotations, il_json, il_annotations = rna3d_motif.load_current_rna3datlas()
    for motif in motifs:
        if motif.loop_type == "hairpin":
            hl_annotations = motif.get_alignments(hl_annotations)
            if motif.include_json:
                motif.get_rna3d_json_sequences(hl_json)
        if motif.loop_type == "internal":
            il_annotations == motif.get_alignments(il_annotations)
            if motif.include_json:
                motif.get_rna3d_json_sequences(il_json)
        motif.get_rna3d_api_sequences()

    rna3d_hairpins = flatten([x.seq_abb_tpls for x in motifs if x.loop_type == "hairpin"])
    rna3d_internals = flatten([x.internal_sequences for x in motifs if x.loop_type == "internal"])

    rna3d_bulges = flatten([x.bulge_sequences for x in motifs if x.loop_type == "internal"])

    # retrieved from tuple created by get_rfam_sequences since they are pre-set
    rfam_hairpins = rfam_sequences_tuples[0]
    # retrieved from tuple created by get_rfam_sequences since they are pre-set
    rfam_internals = rfam_sequences_tuples[1]
    # There are no bulge sequences from Rfam so no need to flatten anything, Rfam bulges file is left empty anyways (it's just there for consistency, if anything changes I'll have to change the rfam_sequences python script)
    rfam_bulges = []

    #################
    # THERE ARE NO EXTRA BULGE SETS TO MAKE SINCE THERE ARE NO BULGES NOTED IN THE RFAM
    # SHOULD THIS EVER CHANGE I WILL NEED TO ADD THEM ABOVE AND CHANGE THE BULGES DUPE_CHECK_WRITE_CSV CALLS
    ################
    both_hairpins = flatten([rna3d_hairpins, rfam_hairpins])
    both_internals= flatten([rna3d_internals, rfam_internals])
    both_bulges = flatten([rna3d_bulges, rfam_bulges])

    write_csv(rna3d_hairpins, "rna3d_hairpins.csv")
    write_csv(rna3d_internals, "rna3d_internals.csv")
    write_csv(rna3d_bulges, "rna3d_bulges.csv")
    write_csv(rfam_hairpins, "rfam_hairpins.csv")
    write_csv(rfam_internals, "rfam_internals.csv")
    write_csv(rfam_bulges, "rfam_bulges.csv")
    write_csv(both_hairpins, "both_hairpins.csv")
    write_csv(both_internals, "both_internals.csv")
    write_csv(both_bulges, "both_bulges.csv")
    update_hexdumbs()


if __name__ == "__main__":
    main()  # main function for getting the sequences from the RNA3DAtlas and the Rfam and writing them to a csv file.
