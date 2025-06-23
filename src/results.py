import sys
from dataclasses import dataclass
import logging
from typing import Literal

logger = logging.getLogger("results")


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

    def __init__(self, id: str, result_list: list[str]) -> None:
        super().__init__(id, result_list)

    @property
    def header(self, classifier: str = "Class") -> str:
        """Returns header string of itself, adapted to currently set algorithm"""
        _header: list[str] = ["ID", classifier, "mfe", "motBracket"]
        return result.separator.join(_header) + "\n"


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
        result_str: str,
        stderr: str,
    ) -> None:
        self.id = name
        self.results = result_str
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

    @results.setter
    def results(self, result_str: str) -> None:
        reslist: list[result_mfe | result_pfc | result_alignment] = []
        split = result_str.strip().split("\n")
        for output in split:
            split_result = output.split("|")
            split_stripped_results = [x.strip() for x in split_result]
            match self.Status:
                case "pfc":
                    res = result_pfc(self.id, split_stripped_results)
                case "mfe":
                    res = result_mfe(self.id, split_stripped_results)
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
        if not initiated:
            if self.stderr.strip():
                logger.warning(self.stderr.strip())
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
