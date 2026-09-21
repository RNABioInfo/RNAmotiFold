import src.RNAmotiFold.results.mfe
import src.RNAmotiFold.results.pfc
import src.RNAmotiFold.results.ali
from typing import Literal, Any
import sys
import logging
from dataclasses import dataclass

logger = logging.getLogger("results")


# List flattening
def flatten(xss: list[list[Any]]) -> list[Any]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]


class algorithm_output:
    """Bigger algorithm output class that mainly holds a list of result objects."""

    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(0, len(v)):
            together.append("{key}={value!r}".format(key=k[i], value=v[i]))
        return f"{classname}({', '.join(together)})"

    def __iter__(self):
        return self

    def __next__(
        self,
    ) -> (
        src.RNAmotiFold.results.mfe.result_mfe
        | src.RNAmotiFold.results.pfc.result_pfc
        | src.RNAmotiFold.results.ali.result_alignment
    ):
        if self._index < len(self.results):
            item = self.results[self._index]
            self._index += 1
            return item
        else:
            self._index = 0
            raise StopIteration

    def __str__(self):
        return "\n".join([self.results[0].header, "\n".join([x.tsv for x in self.results])])

    def __init__(
        self,
        name: str,
        result_str: (
            str
            | list[
                src.RNAmotiFold.results.mfe.result_mfe
                | src.RNAmotiFold.results.pfc.result_pfc
                | src.RNAmotiFold.results.ali.result_alignment
            ]
        ),
        stderr: list[str],
        process_type: Literal["mfe", "pfc", "ali"],
        motif: Literal["hairpin", "internal", "bulge", "all"] = "all",
    ) -> None:
        self.id = name
        self.motif_type = motif
        self.stderr = stderr
        self._index = 0
        self.process_type: Literal["mfe"] | Literal["pfc"] | Literal["ali"] = process_type
        self.results = result_str

    @property
    def results(
        self,
    ) -> list[
        src.RNAmotiFold.results.mfe.result_mfe
        | src.RNAmotiFold.results.pfc.result_pfc
        | src.RNAmotiFold.results.ali.result_alignment
    ]:
        return self._results

    @results.setter
    def results(
        self,
        result: (
            str
            | list[
                src.RNAmotiFold.results.mfe.result_mfe
                | src.RNAmotiFold.results.pfc.result_pfc
                | src.RNAmotiFold.results.ali.result_alignment
            ]
        ),
    ) -> None:
        if isinstance(result, list):
            self._results = result
        else:
            reslist: list[
                src.RNAmotiFold.results.mfe.result_mfe
                | src.RNAmotiFold.results.pfc.result_pfc
                | src.RNAmotiFold.results.ali.result_alignment
            ] = []
            split = result.strip().split("\n")
            match self.process_type:
                case "pfc":
                    pfc_sum = float(sum([float(x.split("|")[1]) for x in split]))
                    for output in split:
                        res = src.RNAmotiFold.results.pfc.result_pfc.from_string(
                            self.id, output, pfc_sum
                        )
                        reslist.append(res)
                case "mfe":
                    for output in split:
                        res = src.RNAmotiFold.results.mfe.result_mfe.from_string(
                            self.id, output
                        )
                        reslist.append(res)
                case "ali":
                    for output in split:
                        res = src.RNAmotiFold.results.ali.result_alignment.from_string(
                            self.id, output
                        )
                        reslist.append(res)
                case _:
                    raise ValueError(f"Invalid result process type detected: {self.process_type}")
            self._results = sorted(reslist)
            if self.process_type == "pfc":
                self.results.reverse()  # pfc has to be flipped because bigger pfc  --> more probable

    @property
    def stderr(self) -> list[str]:
        return self._stderr

    @stderr.setter
    def stderr(self, err: str | list[str]) -> None:
        if isinstance(err, list):
            self._stderr = err
        else:
            if err:
                errlist: list[str] = []
                errlist.append(err.strip())
                self._stderr = errlist
            else:
                self._stderr = []

    # If not initiated function writes a header and then all it's results as csv
    def write_results(self, initiated: bool) -> Literal[True]:
        """Header and results written with this function will be in csv format using the classwide results.separator variable"""
        for err in self.stderr:
            if len(err.strip()) > 0:
                logger.warning(self.id + ": " + err.strip())
        if not initiated:
            sys.stdout.write(self.results[0].header + "\n")
        for result_obj in self.results:
            sys.stdout.write(result_obj.tsv + "\n")
        return True

    @classmethod
    def merge_mfe_outputs(cls, objs: list["algorithm_output"]) -> "algorithm_output":
        """
        Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input
        """
        result_set: set[
            src.RNAmotiFold.results.mfe.result_mfe
            | src.RNAmotiFold.results.pfc.result_pfc
            | src.RNAmotiFold.results.ali.result_alignment
        ] = set()
        for obj in objs:
            for res in obj.results:
                if isinstance(res, src.RNAmotiFold.results.mfe.result_mfe):
                    result_set.add(res)
        sorted_results = sorted(
            list(result_set),
            key=lambda x: (
                x.free_energy
                if isinstance(x, src.RNAmotiFold.results.mfe.result_mfe)
                else 0
            ),
        )
        return cls(
            objs[0].id,
            sorted_results,
            stderr=flatten([x.stderr for x in objs]),
            process_type=objs[0].process_type,
        )


@dataclass
class error:
    id: str
    error: str
