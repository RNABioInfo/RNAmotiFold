import json
from pathlib import Path
import requests
from os import scandir, remove
import json

non_listed_conversions = {"MAD": "A"}

motifs_folder_path = Path(__file__).resolve().parent.joinpath("data", "motifs")
duplicates_json_path = Path(__file__).resolve().parent.joinpath("data", "duplicates.json")


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
            if len(element) > 1:  # if it is a nonstandard base it will have
                element = rna3d_motif.nucleotide_conversion(element)

        return element

    @staticmethod
    def nucleotide_conversion(nucleotide):
        if nucleotide in non_listed_conversions.keys():
            return non_listed_conversions[nucleotide]
        with open("/home/ubuntu/RNAmotiFold/src/data/nucleotide_conversion.json", "r") as file:
            conversion_json = json.load(file)
        return conversion_json[nucleotide]["standard_base"][0]

    @staticmethod
    def reverser(input_list):
        rev_list = []
        for entry in input_list:
            rev_list.append(entry[::-1])
        return rev_list

    @staticmethod
    def get_break(nucleotide_list: list) -> int | None:
        numbers = [rna3d_motif.get_nucleotide_element(x, 4) for x in nucleotide_list]
        chains = [rna3d_motif.get_nucleotide_element(x, 2) for x in nucleotide_list]
        for i in range(len(chains) - 1):
            if chains[i + 1] != chains[i] or abs(int(numbers[i + 1]) - int(numbers[i])) > 5:
                return i
            else:
                pass
        return None

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
    def rna3d_fw(self):
        if hasattr(self, "_rna3d_fw"):
            return self._rna3d_fw
        else:
            return []

    @rna3d_fw.setter
    def rna3d_fw(self, new_data):
        if hasattr(self, "_rna3d_fw"):
            self._rna3d_fw.append(new_data)
        else:
            self._rna3d_fw = []
            self._rna3d_fw.append(new_data)

    @property
    def rna3d_rv(self):
        _bgsu_rv = []
        for element in self.rna3d_fw:
            _bgsu_rv.append(element[::-1])
        return _bgsu_rv

    @property
    def bulge_sequences(self):
        if self.loop_type == "internal":
            return [(x, self.abbreviation) for x in set(self.rna3d_fw) if "$" not in x]
        else:
            raise LookupError(
                "This motif does not have bulge sequences because it is not an internal loop"
            )

    @property
    def internal_sequences(self):
        if self.loop_type == "internal":
            return [(x, self.abbreviation) for x in set(self.rna3d_fw) if "$" in x]
        else:
            raise LookupError(
                "This motif does not have internal sequences because it is not an internal loop"
            )

    @property
    def seq_abb_tpls(self):
        return [(x, self.abbreviation) for x in set(self.rna3d_fw)]

    def load_current_rna3datlas_instances(self):
        match self.loop_type:
            case "hairpin":
                api_response = rna3d_motif.api_call(
                    "http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json"
                )
            case "internal":
                api_response = rna3d_motif.api_call(
                    "http://rna.bgsu.edu/rna3dhub/motifs/release/il/current/json"
                )
        self.motifs = [
            item
            for item in json.loads(api_response.content.decode())
            if item["motif_id"].split(".")[0] in self.rna_3d_atlas_ids
        ]

    def __init__(self, motif_json: dict):
        self.abbreviation = motif_json["abbreviation"]
        self.rna_3d_atlas_ids = motif_json["rna3d_atlas_motif_ids"]
        self.loop_type = motif_json["loop_type"]

    def motifs2csv(self, mot_list):
        return "\n".join([x + f",{self.abbreviation}" for x in set(mot_list)]) + "\n"

    def sequence_from_nucleotide_list(
        self, alignment_nucleotide_list
    ):  # Takes in a list of RNA3DAtlas nucleotide sequences and extracts sequences
        match self.loop_type:
            case "hairpin":
                for alignment in alignment_nucleotide_list:
                    if len(alignment) - 2 > 3:
                        self.rna3d_fw = "".join(
                            rna3d_motif.get_nucleotide_element(x, 3) for x in alignment[1:-1]
                        )
                    else:
                        pass
                        # print(alignment)
            case "internal":
                for alignment in alignment_nucleotide_list:
                    sequence_break = rna3d_motif.get_break(
                        alignment
                    )  # Sequence break is set to none if there wasn't any seq break found or it was below 5
                    if sequence_break is not None:
                        front = alignment[1:sequence_break]
                        back = alignment[sequence_break + 2 : -1]
                        self.rna3d_fw = rna3d_motif.FUSION(
                            "".join(rna3d_motif.get_nucleotide_element(x, 3) for x in front),
                            "".join(rna3d_motif.get_nucleotide_element(x, 3) for x in back),
                        )
                    else:
                        pass
                        # print(alignment)

    def get_rna3d_json_sequences(self):  # Extracts json sequences from the
        for instance in self.motifs:
            alignments = [x for x in instance["alignment"].values()]
            self.sequence_from_nucleotide_list(alignments)

    def get_rna3d_api_sequences(self):
        collection_list = []
        for entry in self.instance_ids:
            api_data = (
                rna3d_motif.api_call(self.rna3datlas_api_call + entry).content.decode().split()
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

    @staticmethod  # collects the rfam_hairpins_fw and rfam_internals_fw
    def get_rfam_sequences():
        rfam_hairpins = motifs_folder_path.joinpath("rfam_hairpins_fw.csv")
        rfam_internals = motifs_folder_path.joinpath("rfam_internals_fw.csv")
        with open(rfam_hairpins) as file:
            hairpins_read = file.readlines()  # type:list[str]
        hairpins = [
            (hairpin.split(",")[0], hairpin.strip().split(",")[1]) for hairpin in hairpins_read
        ]
        with open(rfam_internals) as file2:
            internals_read = file2.readlines()
        internals = [
            (internal.split(",")[0], internal.strip().split(",")[1]) for internal in internals_read
        ]
        return (hairpins, internals)


# Function that takes an input list of tuples(seq,motif abbreviation), checks it for duplicates and writes sequences to csv. Ignores same motif duplicate sequences. Dupes are collected in a .json file.
def dupe_check_write_csv(input_list: list[tuple], filename: str, write_csv: bool = False):
    seen = {}  # type: dict[str,str]
    dupes = {"motif_mode": filename}  # type: dict[str,str|list]
    for tpl in input_list:
        if tpl[0] not in seen.keys():
            seen[tpl[0]] = tpl[1]
        else:
            if seen[tpl[0]] == tpl[1]:
                pass
            else:
                seen[tpl[0]] = seen[tpl[0]].lower()
                if tpl[1] not in dupes.keys():
                    dupes[seen[tpl[0]]] = [tpl[1]]
                else:
                    dupes[seen[tpl[0]]].append(tpl[1])
    if write_csv:
        with open(motifs_folder_path.joinpath(filename), mode="w+") as file:
            file.write("\n".join([x + f",{seen[x]}" for x in seen]))
    return dupes


def flatten(xss):  # Makes a list of lists into a single list with all the items
    return [x for xs in xss for x in xs]


def sequence_reverser(
    seq_abb_tpls: list[tuple[str, str]]
):  # Takes a list of tuples and creates a new list of tuples where the first entry is reversed
    return [(x[0][::-1], x[1]) for x in seq_abb_tpls]


def update_hexdumbs():  # Function iterates of motifs csv files and writes all their hexdums to submodules/RNALoops/Extensions/mot_header.hh
    """main function for updating hexdumbs. Can be used on it's own to update the mot_header.hh file in submodules/RNALoops/Extensions. Does not check for duplicates again, if you want your duplicate.json updated run update_motif_collection.py"""
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
    """Updates motif sequence collection and updates both the duplicates.json file as well as the hexdumb file RNAmotiFold/submodules/RNALoops/Extensions/mot_header.hh"""
    rfam_sequences_tuples = rna3d_motif.get_rfam_sequences()
    motifs = rna3d_motif.load_motif_json()
    for motif in motifs:
        motif.load_current_rna3datlas_instances()
        motif.get_rna3d_json_sequences()
        motif.get_rna3d_api_sequences()

    rna3d_hairpins_fw = flatten([x.seq_abb_tpls for x in motifs if x.loop_type == "hairpin"])
    rna3d_hairpins_rv = sequence_reverser(rna3d_hairpins_fw)
    rna3d_hairpins_both = flatten([rna3d_hairpins_fw, rna3d_hairpins_rv])

    rna3d_internals_fw = flatten(
        [x.internal_sequences for x in motifs if x.loop_type == "internal"]
    )
    rna3d_internals_rv = sequence_reverser(rna3d_internals_fw)
    rna3d_internals_both = flatten([rna3d_internals_fw, rna3d_internals_rv])

    rna3d_bulges_fw = flatten([x.bulge_sequences for x in motifs if x.loop_type == "internal"])
    rna3d_bulges_rv = sequence_reverser(rna3d_bulges_fw)
    rna3d_bulges_both = flatten([rna3d_bulges_fw, rna3d_bulges_rv])

    # retrieved from tuple created by get_rfam_sequences since they are pre-set
    rfam_hairpins_fw = rfam_sequences_tuples[0]
    rfam_hairpins_rv = sequence_reverser(rfam_hairpins_fw)
    rfam_hairpins_both = flatten([rfam_hairpins_fw, rfam_hairpins_rv])

    # retrieved from tuple created by get_rfam_sequences since they are pre-set
    rfam_internals_fw = rfam_sequences_tuples[1]
    rfam_internals_rv = sequence_reverser(rfam_internals_fw)
    rfam_internals_both = flatten([rfam_internals_fw, rfam_internals_rv])

    both_hairpins_fw = flatten([rna3d_hairpins_fw, rfam_hairpins_fw])
    both_hairpins_rv = flatten([rna3d_hairpins_rv, rfam_hairpins_rv])
    both_hairpins_both = flatten([rna3d_hairpins_both, rfam_hairpins_both])

    both_internals_fw = flatten([rna3d_internals_fw, rna3d_internals_fw])
    both_internals_rv = flatten([rna3d_internals_rv, rna3d_internals_rv])
    both_internals_both = flatten([rna3d_internals_both, rfam_internals_both])

    #################
    # THERE ARE NO EXTRA BULGE SETS TO MAKE SINCE THERE ARE NO BULGES NOTED IN THE RFAM
    # SHOULD THIS EVER CHANGE I WILL NEED TO ADD THEM ABOVE AND CHANGE THE BULGES DUPE_CHECK_WRITE_CSV CALLS
    ################
    if Path.is_file(duplicates_json_path):
        remove(duplicates_json_path)
    duplicates = []  # type:list[dict]
    duplicates.append(dupe_check_write_csv(rna3d_hairpins_fw, "rna3d_hairpins_fw.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_hairpins_rv, "rna3d_hairpins_rv.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_hairpins_both, "rna3d_hairpins_both.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_internals_fw, "rna3d_internals_fw.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_internals_rv, "rna3d_internals_rv.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_internals_both, "rna3d_internals_both.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_bulges_fw, "rna3d_bulges_fw.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_bulges_rv, "rna3d_bulges_rv.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_bulges_both, "rna3d_bulges_both.csv"))
    duplicates.append(dupe_check_write_csv(rfam_hairpins_fw, "rfam_hairpins_fw.csv"))
    duplicates.append(dupe_check_write_csv(rfam_hairpins_rv, "rfam_hairpins_rv.csv"))
    duplicates.append(dupe_check_write_csv(rfam_hairpins_both, "rfam_hairpins_both.csv"))
    duplicates.append(dupe_check_write_csv(rfam_internals_fw, "rfam_internals_fw.csv"))
    duplicates.append(dupe_check_write_csv(rfam_internals_rv, "rfam_internals_rv.csv"))
    duplicates.append(dupe_check_write_csv(rfam_internals_both, "rfam_internals_both.csv"))
    duplicates.append(dupe_check_write_csv(both_hairpins_fw, "both_hairpins_fw.csv"))
    duplicates.append(dupe_check_write_csv(both_hairpins_rv, "both_hairpins_rv.csv"))
    duplicates.append(dupe_check_write_csv(both_hairpins_both, "both_hairpins_both.csv"))
    duplicates.append(dupe_check_write_csv(both_internals_fw, "both_internals_fw.csv"))
    duplicates.append(dupe_check_write_csv(both_internals_rv, "both_internals_rv.csv"))
    duplicates.append(dupe_check_write_csv(both_internals_both, "both_internals_both.csv"))

    # THERE ARE ACTUALLY NO BULGES IN THE RFAM SO IF I EVER ADD ONE OF THEIR MOTIFS THAT INCLUDES A BULGE OR THEY UPDATE THEIR DATABASE ADD THEM BELOW HERE
    duplicates.append(dupe_check_write_csv(rna3d_bulges_fw, "both_bulges_fw.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_bulges_rv, "both_bulges_rv.csv"))
    duplicates.append(dupe_check_write_csv(rna3d_bulges_both, "both_bulges_both.csv"))

    with open(duplicates_json_path, "w+") as file:
        file.write(json.dumps(duplicates))

    update_hexdumbs()


if __name__ == "__main__":
    main()  # main function for getting the sequences from the RNA3DAtlas and the Rfam and writing them to a csv file.
