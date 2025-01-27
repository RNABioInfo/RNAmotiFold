import requests
from Bio import AlignIO
from io import StringIO
from pathlib import Path


class rfam_motif:
    def __init__(self, abb: str, id: str, lower_bound=None, upper_bound=None):
        self.abbreviation = abb
        self.id = id
        self.lower = lower_bound
        self.upper = upper_bound

    def get_rfam_alignments(self):
        i = 0
        while i < 6:
            answer = requests.get(
                f"https://rfam.org/motif/{self.id}/alignment?acc={self.id}&format=stockholm&download=0"
            )
            if answer.status_code == 200:
                self.alignment = list(
                    AlignIO.parse(StringIO(answer.content.decode()), format="stockholm")
                )[0]
                break
            else:
                i += 1
        if i == 6:
            raise ConnectionError(f"Could not retrieve rfam alignment for {self.motif_name}.")

    def extract_rfam_sequences(self):
        for entry in self.alignment._records:  # type:ignore
            seq = [x.upper() for x in list(entry.seq[self.lower : self.upper]) if x != "-"]
            if "N" not in seq:
                self.seqs = "".join(seq)
            else:
                pass

    @property
    def seqs(self):
        return list(set(self._seqs))

    @seqs.setter
    def seqs(self, newseq: str):
        if hasattr(self, "seqs"):
            self._seqs.append(newseq)
        else:
            self._seqs = [newseq]

    def motifs2csv(self, mot_list):
        return "\n".join([x + f",{self.abbreviation}" for x in mot_list]) + "\n"

    def motifs_rev2csv(self, mot_list):
        return "\n".join([x[::-1] + f",{self.abbreviation}" for x in mot_list]) + "\n"


def main():
    T = rfam_motif(abb="T", id="RM00024", lower_bound=21, upper_bound=30)
    G = rfam_motif(abb="G", id="RM00008", lower_bound=35, upper_bound=39)
    U = rfam_motif(abb="U", id="RM00029", lower_bound=24, upper_bound=28)
    mots = [T, G, U]  # type:list[rfam_motif]
    initialized = False
    for mot in mots:
        mot.get_rfam_alignments()
        mot.extract_rfam_sequences()
        if not initialized:
            write_mode = "w+"
            initialized = True
        else:
            write_mode = "a"
        with open(
            f"{Path(__file__).resolve().parent.joinpath("data","motifs","rfam_hairpins_fw.csv")}",
            write_mode,
        ) as file:
            file.write(mot.motifs2csv(mot.seqs))


if __name__ == "__main__":
    main()
