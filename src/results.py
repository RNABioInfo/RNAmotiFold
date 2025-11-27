import sys
from dataclasses import dataclass
import logging
from typing import Literal,Any
import re

logger = logging.getLogger("results")

#List flattening
def flatten(xss:list[list[Any]]) -> list[Any]:
    '''
    Used to flatten a lists of lists into a single list
    '''
    return [x for xs in xss for x in xs]


class result:
    separator: str = ","

    def __init__(self, id: str, result_list: list[str]) -> None:
        self.id: str = id
        self.cols: list[str] = result_list

    def __str__(self) -> str:
        return self.tsv
    
    @property
    def tsv(self) -> str:
        """Returns tsv string of itself"""
        return self.id + result.separator + result.separator.join(map(str, self.cols)) + "\n"


class result_mfe(result):

    def __eq__(self,other:'result_mfe') -> bool: #type:ignore
        return self.free_energy == other.free_energy and self.mot_bracket == other.mot_bracket and self.classifier == other.classifier

    def __hash__(self) -> int:
        return hash((self.free_energy,self.mot_bracket,self.classifier))
    
    @classmethod
    def merge_structures(cls,compatibles:list['result_mfe']) ->'result_mfe|None':
        compatibles.sort(key= lambda x: x.classifier[0]) #sort list in place alternative would be new = sorted(compatibles,key=...)
        base_structure = list(compatibles[0].dot_bracket)
        insertions:set[int] = set()
        motifs:set[tuple[str,str]] = set()
        for result in compatibles:
                motif= result.classifier[0]
                locations = list(result_mfe.find_all(result.mot_bracket,motif))
                for loc in locations:
                    if loc in insertions and base_structure[loc] != motif:
                        motifs.add((motif.lower(),result.motif_type))
                        logger.warning(f"Overlap detected during merge at position {loc} between motifs {base_structure[loc]} and {motif} on sequence: {result.id}, marking as lowercase and inserting {base_structure[loc].lower()}")
                        base_structure[loc] = base_structure[loc].lower()
                    else:
                        motifs.add((motif,result.motif_type))
                        base_structure[loc] = motif
                    insertions.add(loc)
        merged_bracket = "".join(base_structure)
        foundslist:list[tuple[int,tuple[str,str]]] = []
        for m in motifs:
            founds = re.finditer(f"[()]{m[0]}+",merged_bracket)
            for f in founds:
                foundslist.append((f.start(),m))
        foundslist.sort(key=lambda tup:tup[0])
        new_classifier = result_mfe.build_new_classifier([x[1] for x in foundslist])
        if merged_bracket not in [x.mot_bracket for x in compatibles]:
            return cls(compatibles[0].id+"_merged",[new_classifier,str(compatibles[0].free_energy),merged_bracket],"all")
        else:
            return None

    @staticmethod
    def build_new_classifier(foundslist:list[tuple[str,str]]) -> str:
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
    def find_all(a_str:str, sub:str):
        start = 0
        while True:
            start = a_str.find(sub, start)
            if start == -1: return
            yield start
            start += 1 # use start += 1 to find overlapping matches

    @staticmethod
    def get_compatible_structures(struc_list:list['result_mfe']) -> list[list[int]]:
        collecting:dict[int,list[int]] = {}
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
        compatible = [[k]+v for k,v in collecting.items()]
        return compatible

    def __init__(self, id: str, result_list: list[str],motif_type:Literal["hairpin","internal","bulge","all"]) -> None:
        super().__init__(id, result_list)
        self.classifier = self.cols[0]
        self.free_energy = self.cols[1]
        self.mot_bracket = self.cols[2]
        self.dot_bracket = self.cols[2]
        self.motif_type:Literal["hairpin","internal","bulge","all"] = motif_type
        
    @property
    def dot_bracket(self) -> str:
        return self._dot_bracket
    
    @dot_bracket.setter
    def dot_bracket(self,structure_string:str) -> None:
        for c in set(self.classifier):
            structure_string = structure_string.replace(c,".")
        self._dot_bracket = structure_string

    @property
    def header(self, classifier: str = "Class") -> str:
        """Returns header string of itself, adapted to currently set algorithm"""
        _header: list[str] = ["ID", classifier, "mfe", "motBracket"]
        return result.separator.join(_header) + "\n"

    @property
    def classifier(self) -> str:
        return self._classifier
    
    @classifier.setter
    def classifier(self,class_string:str):
        self._classifier = class_string

    @property
    def free_energy(self) -> int:
        return self._free_energy
    
    @free_energy.setter
    def free_energy(self, new_energy:str) -> None:
        self._free_energy = int(new_energy)

    @property
    def mot_bracket(self) -> str:
        return self._mot_bracket
    
    @mot_bracket.setter
    def mot_bracket(self,motbracket_string:str) -> None:
        self._mot_bracket = motbracket_string


class result_pfc(result):
    def __init__(self, id: str, result_list: list[str]) -> None:
        super().__init__(id, result_list)

    @property
    def header(self, classifier: str = "Class") -> str:
        """Returns header string of itself, adapted to currently set algorithm"""
        _header: list[str] = ["ID", classifier, "pfc", "Probability"]
        return result.separator.join(_header) + "\n"


class result_alignment(result):
    """Dummy class for compatibility, fill out later"""

    @property
    def header(self) -> str:
        return "bruh"


@dataclass
class error:
    id: str
    error: str


class algorithm_output:
    """Bigger algorithm output class that mainly holds a list of result objects."""

    # Result type
    _Status: Literal["pfc", "mfe", "alignment"]

    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(0, len(v)):
            together.append("{key}={value!r}".format(key=k[i], value=v[i]))
        return f"{classname}({', '.join(together)})"

    def __iter__(self):
        return self

    def __next__(self) -> result:
        if self._index < len(self.results):
            item = self.results[self._index]
            self._index += 1
            return item
        else:
            self._index = 0
            raise StopIteration

    def __str__(self):
        return self.results[0].header + "".join([str(x) for x in self.results])

    def __init__(
        self,
        name: str,
        result_str: str|list[result_mfe|result_pfc|result_alignment],
        stderr: list[str],
        motif_type:Literal["hairpin","internal","bulge","all"] = "all",
    ) -> None:
        self.id = name
        self.results = (result_str,motif_type)
        self.stderr = stderr
        self._index = 0

    # Formats results from the mgapc output formatting to a list of result objects
    @property
    def Status(self) -> str:
        return self._Status

    @Status.setter
    def Status(self, status: Literal["pfc", "mfe", "alignment"]):
        self._Status = status

    @property
    def results(self) -> list[result_mfe | result_pfc | result_alignment]:
        return self._results

    @property
    def stderr(self) -> list[str]:
        return self._stderr
    
    @stderr.setter
    def stderr(self,err:str|list[str]) -> None:
        if isinstance(err,list):
            self._stderr = err
        else:
            if err:
                errlist:list[str] = []
                errlist.append(err.strip())
                self._stderr = errlist
            else:
                self._stderr = []

    @results.setter
    def results(self, result:tuple[str|list[result_mfe|result_pfc|result_alignment],Literal["hairpin","internal","bulge","all"]]) -> None:
        if isinstance(result[0],list):
            self._results = result[0]
            return None
        reslist: list[result_mfe | result_pfc | result_alignment] = []
        split = result[0].strip().split("\n")
        for output in split:
            split_result = output.split("|")
            split_stripped_results = [x.strip() for x in split_result]
            match self.Status:
                case "pfc":
                    res = result_pfc(self.id, split_stripped_results)
                case "mfe":
                    res = result_mfe(self.id, split_stripped_results,result[1])
                case "alignment":
                    res = result_alignment(self.id, split_stripped_results)
                case _:
                    raise ValueError(f"Invalid result status detected: {self.Status}")
            reslist.append(res)
        self._results = reslist
        if self.Status == "pfc":
            self.calculate_pfc_probabilities()

    # If not initiated function writes a header and then all it's results as csv
    def write_results(self, initiated: bool) -> Literal[True]:
        """Header and results written with this function will be in csv format using the classwide results.separator variable"""
        for err in self.stderr:
            if len(err.strip()) > 0:
                logger.warning(self.id+": "+err.strip())
        if not initiated:
            sys.stdout.write(self.results[0].header)
        for result_obj in self.results:
            sys.stdout.write(result_obj.tsv)
        return True

    # Calculates partition function probabilities for me so I don't have to manually do it every time for all the outputs
    def calculate_pfc_probabilities(self) -> None:
        pfc_list: list[float] = []
        for result_obj in self.results:
            pfc_val = float(result_obj.cols[-1])
            pfc_list.append(pfc_val)
        pfc_sum = sum(pfc_list)
        for result_obj in self.results:
            result_obj.cols.append(str(round(float(result_obj.cols[-1]) / pfc_sum, 4)))
    
    @classmethod
    def merge_outputs(cls,objs:list['algorithm_output']) -> 'algorithm_output':
        '''
        Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input
        '''
        result_set:set[result_mfe |result_pfc | result_alignment] = set()
        for obj in objs:
            for res in obj.results:
                if isinstance(res,result_mfe):
                    result_set.add(res)
        return cls(objs[0].id,list(result_set),stderr=flatten([x.stderr for x in objs]))