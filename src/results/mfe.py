from results import result
import re


class result_mfe(result):
    """Subclass of result for mfe results, has extra attributes for free energy and motBracket structure. Also implements comparison and hashing for mering structures in single motif mode"""

    def __init__(self, id: str, classifier: str, free_energy: str, mot_bracket: str) -> None:
        super().__init__(id, classifier)
        self.free_energy = free_energy
        self.motBracket = mot_bracket #This variable gets special treatment so our outputs looks nice

    # Special Dunder Method for hashing and comparing mfe results, used for merging structure in single motif mode
    def __eq__(self, other: object) -> bool:
        if isinstance(other, result_mfe):
            return (
                self.free_energy == other.free_energy
                and self.motBracket == other.motBracket
                and self.classifier == other.classifier
            )
        else:
            raise NotImplementedError(f"Cannot compare result_mfe with {type(other)}")

    def __ne__(self, other: object) -> bool:
        if isinstance(other, result_mfe):
            return not self.__eq__(other)
        else:
            raise NotImplementedError(f"Cannot compare result_mfe with {type(other)}")

    def __lt__(self, other: "result_mfe") -> bool:
        return self.free_energy < other.free_energy

    def __le__(self, other: "result_mfe") -> bool:
        return self.free_energy <= other.free_energy

    def __gt__(self, other: "result_mfe") -> bool:
        return self.free_energy > other.free_energy

    def __ge__(self, other: "result_mfe") -> bool:
        return self.free_energy >= other.free_energy

    def __hash__(self) -> int:
        return hash((self.free_energy, self.motBracket, self.classifier))

    @property
    def dot_bracket(self) -> str:
        return self._dot_bracket

    @dot_bracket.setter
    def dot_bracket(self, structure_string: str) -> None:
        for c in set(self.classifier):
            structure_string = structure_string.replace(c, ".")
        self._dot_bracket = structure_string

    @property
    def free_energy(self) -> float:
        return self._free_energy  # type: ignore Can't be anything but float because setter only permits float

    @free_energy.setter
    def free_energy(self, new_energy: str) -> None:
        self._free_energy = float(int(new_energy) / 100)

    @classmethod
    def from_string(cls, id: str, result_string: str) -> "result_mfe":
        split_result = result_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(
            id=id,
            classifier=split_stripped_results[0],
            free_energy=split_stripped_results[1],
            mot_bracket=split_stripped_results[2],
        )

    # Implementation of structure merging for single motif mode
    @classmethod
    def merge_structures(cls, compatibles: list["result_mfe"]) -> "result_mfe|None":
        compatibles.sort(
            key=lambda x: x.classifier[0]
        )  # sort list in place alternative would be new = sorted(compatibles,key=...)
        base_structure = list(compatibles[0].dot_bracket)
        insertions: set[int] = set()
        motifs: set[tuple[str, str]] = set()
        for result in compatibles:
            motif = result.classifier[0]
            locations = list(result_mfe.find_all(result.motBracket, motif))
            for loc in locations:
                if loc in insertions and base_structure[loc] != motif:
                    motifs.add((motif.lower(), result.motif_type))
                    base_structure[loc] = base_structure[loc].lower()
                else:
                    motifs.add((motif, result.motif_type))
                    base_structure[loc] = motif
                insertions.add(loc)
        merged_bracket = "".join(base_structure)
        foundslist: list[tuple[int, tuple[str, str]]] = []
        for m in motifs:
            founds = re.finditer(f"[()]{m[0]}+", merged_bracket)
            for f in founds:
                foundslist.append((f.start(), m))
        foundslist.sort(key=lambda tup: tup[0])
        new_classifier = result_mfe.build_new_classifier([x[1] for x in foundslist])
        if merged_bracket not in [x.motBracket for x in compatibles]:
            return cls(
                id=compatibles[0].id + "_merged",
                classifier=new_classifier,
                free_energy=str(compatibles[0].free_energy),
                mot_bracket=merged_bracket,
            )
        else:
            return None

    @staticmethod
    def build_new_classifier(foundslist: list[tuple[str, str]]) -> str:
        new_classifier = ""
        for found in foundslist:
            if found[1] == "hairpin":
                new_classifier += found[0]
            elif found[1] == "bulge":
                new_classifier += found[0]
            elif found[1] == "internal":
                new_classifier += found[0]
                foundslist.reverse()
                foundslist.remove(found)
                foundslist.reverse()
            else:
                raise ValueError("Invalid motif type detected during merge")
        return new_classifier

    @staticmethod
    def find_all(a_str: str, sub: str):
        start = 0
        while True:
            start = a_str.find(sub, start)
            if start == -1:
                return
            yield start
            start += 1  # use start += 1 to find overlapping matches

    @staticmethod
    def get_compatible_structures(struc_list: list["result_mfe"]) -> list[list[int]]:
        collecting: dict[int, list[int]] = {}
        for i in range(len(struc_list)):
            for j in range(len(struc_list)):
                if i >= j:
                    continue
                first = struc_list[i]
                second = struc_list[j]
                if first.dot_bracket == second.dot_bracket:
                    if first not in collecting.keys():
                        collecting[i] = [j]
                    else:
                        collecting[i].append(j)
        compatible = [[k] + v for k, v in collecting.items()]
        return compatible
