import sys
from dataclasses import dataclass
import logging
from typing import Literal,Any, NamedTuple
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

    def __init__(self,id:str,classifier:str,motif_type:Literal["hairpin","internal","bulge","all"]) -> None:
        self.values:dict[str,str|int|float] = dict()
        self.values["ID"] = id
        self.values["class"] = classifier
        self.motif_type:Literal["hairpin","internal","bulge","all"] = motif_type

    @property
    def id(self) -> str:
        return self.values["ID"] #type:ignore Can't be anything but string because setter only permits string
    
    @id.setter
    def id(self, new_id:str) -> None:
        self.values["ID"] = new_id
        
    @property
    def classifier(self) -> str:
        return self.values["class"] #type:ignore Can't be anything but string because setter only permits string
    
    @classifier.setter
    def classifier(self, new_classifier:str) -> None:
        self.values["class"] = new_classifier

    def __str__(self) -> str:
        return self.tsv
    
    def add_coloumn(self,name:str,value:str|int|float) -> None:
        self.values[name] = value

    @property
    def tsv(self) -> str:
        """Returns tsv string of itself"""
        return result.separator.join([str(x) for x in self.values.values()]) + "\n"
    
    @property
    def header(self) -> str:
        """Returns header string of itself, adapted to currently set algorithm"""
        return result.separator.join(self.values.keys()) + "\n"


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
            return cls(id=compatibles[0].id+"_merged",classifier=new_classifier,free_energy=str(compatibles[0].free_energy),mot_bracket=merged_bracket,motif_type=compatibles[0].motif_type)
        else:
            return None

    @classmethod
    def from_string(cls,id:str,result_string:str,motif_type:Literal["hairpin","internal","bulge","all"]) -> 'result_mfe':
        split_result = result_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(id=id,classifier=split_stripped_results[0],free_energy=split_stripped_results[1],mot_bracket=split_stripped_results[2],motif_type=motif_type)

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

    def __init__(self, id: str, classifier: str, free_energy: str, mot_bracket: str, motif_type:Literal["hairpin","internal","bulge","all"]) -> None:      
        super().__init__(id,classifier,motif_type)
        self.free_energy = free_energy
        self.mot_bracket = mot_bracket

    @property
    def dot_bracket(self) -> str:
        return self._dot_bracket
    
    @dot_bracket.setter
    def dot_bracket(self,structure_string:str) -> None:
        for c in set(self.classifier):
            structure_string = structure_string.replace(c,".")
        self._dot_bracket = structure_string

    @property
    def mot_bracket(self) -> str:
        return self.values["motBracket"] #type:ignore Can't be anything but string because setter only permits string

    @mot_bracket.setter
    def mot_bracket(self,motbracket_string:str) -> None:
        self.values["motBracket"] = motbracket_string

    @property
    def free_energy(self) -> float:
        return self.values["mfe"] #type:ignore Can't be anything but float because setter only permits float

    @free_energy.setter
    def free_energy(self, new_energy:str) -> None:
        self.values["mfe"] = float(int(new_energy)/100)


class result_pfc(result):
    def __init__(self, id: str, classifier:str, pfc_value:str|int|float,motif_type:Literal["hairpin","internal","bulge","all"]) -> None:
        super().__init__(id, classifier,motif_type)
        self.pfc_value = pfc_value
    
    @property
    def pfc_value(self) -> float:
        return float(self.values["pfc"])
    
    @pfc_value.setter
    def pfc_value(self, new_pfc:str|float|int) -> None:
        self.values["pfc"] = new_pfc

    @property
    def probability(self) -> float:
        return self.values["probability"] #type:ignore Can't be anything but float because setter only permits float
    
    @probability.setter
    def probability(self, prob:float) -> None:
        self.values["probability"] = prob

    @classmethod
    def from_string(cls,id:str,result_string:str,motif_type:Literal["hairpin","internal","bulge","all"]) -> 'result_pfc':
        split_result = result_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(id=id,classifier=split_stripped_results[0],pfc_value=split_stripped_results[1],motif_type=motif_type)

    def set_probability(self,pfc_sum:float) -> None:
        self.probability = round(self.pfc_value / pfc_sum, 4)

class alignment_score(NamedTuple):
    energy: float
    covariance: float
    motif: float
    overall:float

    @classmethod
    def from_string(cls,score_string:str):
        nums = [float(x) for x in re.findall(r'-?\d+(?:\.\d+)?', score_string)]
        return cls(overall=nums[0],energy=nums[1],covariance=nums[2],motif=nums[3])
    
    def to_string(self,sep:str):
        return sep.join([str(self.overall),str(self.energy),str(self.covariance),str(self.motif)])

class result_alignment:
    """Dummy class for compatibility, fill out later"""

    def __init__(self,id:str,score:alignment_score,motBracket:str):
        self.id= id
        self.score:alignment_score = score
        self.motBracket:str = motBracket

    @property
    def overall_score(self)->float:
        return self.score.overall
    
    @property
    def energy(self)->float:
        return self.score.energy
    
    @property
    def covariance(self)->float:
        return self.score.covariance
    
    @property
    def motif_score(self) -> float:
        return self.score.motif

    @property
    def header(self) -> str:
       return result.separator.join(["ID","Total Score","Free Energy","Covariance Score","Motif Score","MotBracket"]) + "\n"
    
    @classmethod
    def from_string(cls,id:str,results_string:str) -> 'result_alignment':
        split_result = results_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(id=id, score=alignment_score.from_string(split_stripped_results[0]), motBracket=split_stripped_results[1])

    @property
    def tsv(self):
        return result.separator.join([self.id,str(self.score.overall),str(self.score.energy),str(self.score.covariance),str(self.score.motif),self.motBracket]) +"\n"

    def add_coloumn(self,name:str,value:str|int|float) -> None:
        self.name = value


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

    def __next__(self) -> result|result_alignment:
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
        motif:Literal["hairpin","internal","bulge","all"] = "all",
    ) -> None:
        self.id = name
        self.motif_type = motif
        self.results = result_str
        self.stderr = stderr
        self._index = 0

    @property
    def motif_type(self) -> Literal["hairpin","internal","bulge","all"]:
        return self._motif_type
    
    @motif_type.setter
    def motif_type(self,new_type:Literal["hairpin","internal","bulge","all"]) -> None:
        self._motif_type:Literal["hairpin","internal","bulge","all"] = new_type

    # Formats results from the mgapc output formatting to a list of result objects
    @property
    def Status(self) -> str:
        return self._Status

    @Status.setter
    def Status(self, status: Literal["pfc", "mfe", "alignment"]) -> None:
        self._Status = status

    @property
    def results(self) -> list[result_mfe | result_pfc | result_alignment]:
        return self._results

    @results.setter
    def results(self, result:str|list[result_mfe|result_pfc|result_alignment]) -> None:
        if isinstance(result,list):
            self._results = result
        else:
            reslist: list[result_mfe | result_pfc | result_alignment] = []
            split = result.strip().split("\n")
            match self.Status:
                case "pfc":
                    for output in split:
                        res = result_pfc.from_string(self.id,output,self.motif_type)
                        reslist.append(res)
                case "mfe":
                    for output in split:
                        res = result_mfe.from_string(self.id,output,self.motif_type)
                        reslist.append(res)
                case "alignment":
                    for output in split:
                        res = result_alignment.from_string(self.id,output)
                        reslist.append(res)
                case _:
                    raise ValueError(f"Invalid result status detected: {self.Status}")         
            self._results = sorted(reslist, key=lambda x: x.free_energy if isinstance(x,result_mfe) else x.pfc_value if isinstance(x,result_pfc) else x.overall_score if isinstance(x,result_alignment) else len(x.id))
            if self.Status == "pfc":
                self.results.reverse() #pfc has to be flipped because bigger pfc  --> more probable
                self.add_pfc_probabilities()

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

    def add_column_to_all(self,name:str,value:str|int|float) -> None:
        for res in self.results:
            res.add_coloumn(name,value)

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
    def add_pfc_probabilities(self) -> None:
        if any(not isinstance(x,result_pfc) for x in self.results):
            raise ValueError("Unable to add probabilities, not all results are of type pfc")
        else:
            pfc_sum:float = float(sum([result_obj.pfc_value for result_obj in self.results])) #type:ignore 
            for result_obj in self.results:
                result_obj.set_probability(pfc_sum) #type:ignore Always pfc cause we check beforehand

    @classmethod
    def merge_mfe_outputs(cls,objs:list['algorithm_output']) -> 'algorithm_output':
        '''
        Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input
        '''
        result_set:set[result_mfe |result_pfc|result_alignment] = set()
        for obj in objs:
            for res in obj.results:
                if isinstance(res,result_mfe):
                    result_set.add(res)
        sorted_results = sorted(list(result_set),key=lambda x: x.free_energy if isinstance(x,result_mfe) else 0)
        return cls(objs[0].id,sorted_results,stderr=flatten([x.stderr for x in objs]))