import sys
from dataclasses import dataclass


class result:

    _separator: str = ","
    _algorithm: str = "motmfepretty"

    def __str__(self):
        return self.tsv

    def __init__(self, id: str, result_list: list) -> None:
        self.id = id
        self.cols = result_list

    @property
    def tsv(self):
        """Returns tsv string of itself"""
        return self.id + result._separator + result._separator.join(map(str, self.cols)) + "\n"

    @property
    def separator(self):
        return self._separator

    @property
    def header(self):
        """Returns header string of itself, adapted to currently set algorithm"""
        match self._algorithm:
            case "motmfepretty" | "motmfepretty_subopt":
                self._header = ["ID", "motifs", "mfe", "motBracket"]
            case (
                "mothishape_h"
                | "mothishape_m"
                | "mothishape_b"
                | "mothishape_h_subopt"
                | "mothishape_m_subopt"
                | "mothishape_b_subopt"
            ):
                self._header = ["ID", "motiCes", "mfe", "motBracket"]
            case "motshapeX" | "motshapeX_subopt":
                self._header = ["ID", "mosh", "mfe", "motBracket"]
            case "mothishape_h_pfc" | "mothishape_b_pfc" | "mothishape_m_pfc":
                self._header = ["ID", "MotiCes", "pfc", "probability"]
            case "motshapeX_pfc":
                self._header = ["ID", "mosh", "pfc", "probability"]
            case "motpfc":
                self._header = ["ID", "motifs", "pfc", "probability"]
        return result._separator.join(self._header) + "\n"


@dataclass
class error:
    id: str
    error: str


class algorithm_output:
    """Class specifically made for reading mgapc outputs."""

    # Class wide variables that will be the same for all algorithm outputs at the same time.
    # These get set with the _calibrate_results function during auto_run() in the main bgap_rna class
    pfc = False
    filtering = False

    def __iter__(self):
        return self

    def __next__(self) -> result:
        if self._index < len(self.results):
            item = self.results[self._index]
            self._index += 1
            return item
        else:
            raise StopIteration

    def __str__(self):
        return self.results[0].header + "".join([str(x) for x in self.results])

    def __init__(self, name: str, result_str: str, time_str: str):
        self.id = name
        self.results = self._format_results(result_str)  # type:list[result]
        self.time_str = time_str
        self._index = 0
        if algorithm_output.pfc:
            self.calculate_pfc_probabilities()
            if algorithm_output.filtering:
                self._filter_out_low_probability()

    @property
    def time_str(self):
        return self._time_str

    @time_str.setter
    def time_str(self, time: str):
        if time == "":
            self._time_str = None
        else:
            self._time_str = time

    # Formats results from the mgapc output formatting to a list of result objects
    def _format_results(self, result_str: str) -> list[result]:
        reslist = []
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
            sys.stdout.write(self.results[0].header)
        for result_obj in self.results:
            sys.stdout.write(result_obj.tsv)
        if self.time_str is not None:
            sys.stderr.write(f"{self.id}: {self.time_str}")
        return True

    # Calculates partition function probabilities for me so I don't have to manually do it every time for all the outputs
    def calculate_pfc_probabilities(self) -> None:
        pfc_list = []
        for result_obj in self.results:
            pfc_val = float(result_obj.cols[-1])
            pfc_list.append(pfc_val)
        pfc_sum = sum(pfc_list)
        for result_obj in self.results:
            result_obj.cols.append(round(float(result_obj.cols[-1]) / pfc_sum, 4))

    # Filters out partition function results that have probabilities below 0.0001
    def _filter_out_low_probability(self):
        """Filter out low probability results. Probabilities are calculated with these included"""
        self.results = [result for result in self.results if result.cols[-1] > 0.0001]
