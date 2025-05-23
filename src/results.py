import sys
from dataclasses import dataclass
import logging

logger = logging.getLogger("results")


class result:

    separator: str = ","
    algorithm: str = "RNAmotiFold"

    def __str__(self) ->str:
        return self.tsv

    def __init__(self, id: str, result_list: list[str]) -> None:
        self.id = id
        self.cols = result_list

    @property
    def tsv(self)->str:
        """Returns tsv string of itself"""
        return self.id + result.separator + result.separator.join(map(str, self.cols)) + "\n"

    @property
    def header(self):
        """Returns header string of itself, adapted to currently set algorithm"""
        match self.algorithm:
            case "RNAmotiFold" | "RNAmotiFold_subopt":
                _header = ["ID", "motifs", "mfe", "motBracket"]
            case "RNAmotiCes" | "RNAmotiCes_subopt":
                _header = ["ID", "motiCes", "mfe", "motBracket"]
            case "RNAmoSh" | "RNAmoSh_subopt":
                _header = ["ID", "moSh", "mfe", "motBracket"]
            case "RNAmotiCes_pfc":
                _header = ["ID", "motiCes", "pfc", "probability"]
            case "RNAmoSh_pfc":
                _header = ["ID", "moSh", "pfc", "probability"]
            case "RNAmotiFold_pfc":
                _header = ["ID", "motifs", "pfc", "probability"]
            case _:
                _header = [f"col{self.cols.index(x)}" for x in self.cols]
                _header.insert(0, "ID")
        return result.separator.join(_header) + "\n"





@dataclass
class error:
    id: str
    error: str


class algorithm_output:
    """Bigger algorithm output class that mainly holds a list of result objects."""

    # Class wide variables that will be the same for all algorithm outputs at the same time.
    # These get set with the _calibrate_results function during auto_run() in the main bgap_rna class
    pfc = False

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

    def __init__(self, name: str, result_str: str, stderr: str):
        self.id = name
        self.results: list[result] = self._format_results(result_str)
        self.stderr = stderr
        self._index = 0
        if algorithm_output.pfc:
            self.calculate_pfc_probabilities()

    # Formats results from the mgapc output formatting to a list of result objects
    def _format_results(self, result_str: str) -> list[result]:
        reslist:list[result] = []
        split = result_str.strip().split("\n")
        for output in split:
            split_result = output.split("|")
            split_stripped_results = [x.strip() for x in split_result]
            res = result(self.id, split_stripped_results)
            reslist.append(res)
        return reslist

    # If not initiated function writes a header and then all it's results as csv
    def write_results(self, initiated: bool):
        """Header and results written with this function will be in csv format using the classwide results.separator variable"""
        if not initiated:
            if self.stderr.strip():
                logger.warning(self.stderr.strip())
            sys.stdout.write(self.results[0].header)
        for result_obj in self.results:
            sys.stdout.write(result_obj.tsv)
        return True

    # Calculates partition function probabilities for me so I don't have to manually do it every time for all the outputs
    def calculate_pfc_probabilities(self) -> None:
        pfc_list:list[float] = []
        for result_obj in self.results:
            pfc_val = float(result_obj.cols[-1])
            pfc_list.append(pfc_val)
        pfc_sum = sum(pfc_list)
        for result_obj in self.results:
            result_obj.cols.append(str(round(float(result_obj.cols[-1]) / pfc_sum, 4)))
