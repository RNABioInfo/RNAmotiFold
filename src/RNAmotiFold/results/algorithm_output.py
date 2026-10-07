from RNAmotiFold.results.mfe import ResultMFE
from RNAmotiFold.results.pfc import ResultPFC
from RNAmotiFold.results.ali import ResultAlignment
from typing import Literal, Any
import sys
import logging
from dataclasses import dataclass

logger = logging.getLogger(__name__)


# List flattening
def _flatten(xss: list[list[Any]]) -> list[Any]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]


class AlgorithmOutput:
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
    ) -> ResultMFE | ResultPFC | ResultAlignment:
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
        result_str: str | list[ResultMFE] | list[ResultPFC] | list[ResultAlignment],
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
    ) -> list[ResultMFE] | list[ResultPFC] | list[ResultAlignment]:
        return self._results

    @results.setter
    def results(
        self,
        result: str | list[ResultMFE] | list[ResultPFC] | list[ResultAlignment],
    ) -> None:
        if isinstance(result, list):
            self._results = result
        else:
            split = result.strip().split("\n")
            match self.process_type:
                case "pfc":
                    reslist_pfc: list[ResultPFC] = []
                    pfc_sum = float(sum([float(x.split("|")[1]) for x in split]))
                    for output in split:
                        res_pfc: ResultPFC = ResultPFC._from_string(
                            id=self.id, result_string=output, pfc_sum=pfc_sum
                        )
                        reslist_pfc.append(res_pfc)
                    self._results = sorted(reslist_pfc)
                    self.results.reverse()  # pfc has to be flipped because bigger pfc  --> more probable

                case "mfe":
                    reslist_mfe: list[ResultMFE] = []
                    for output in split:
                        res_mfe: ResultMFE = ResultMFE._from_string(id=self.id, result_string=output)
                        reslist_mfe.append(res_mfe)
                    self._results = sorted(reslist_mfe)

                case "ali":
                    reslist_ali: list[ResultAlignment] = []
                    for output in split:
                        res_ali: ResultAlignment = ResultAlignment._from_string(
                            id=self.id, results_string=output
                        )
                        reslist_ali.append(res_ali)
                    self._results = sorted(reslist_ali)

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
            logger.debug("Starting result writing")
            sys.stdout.write(self.results[0].header + "\n")
        for result_obj in self.results:
            sys.stdout.write(result_obj.tsv + "\n")
        return True

    @classmethod
    def _merge_mfe_outputs(cls, objs: list["AlgorithmOutput"]) -> "AlgorithmOutput":
        """
        Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input
        """
        logger.debug(f"Merging results for {objs[0].id}")
        result_set: set[ResultMFE] = set()
        for obj in objs:
            for res in obj.results:
                if isinstance(res, ResultMFE):
                    result_set.add(res)
        sorted_results = sorted(list(result_set), key=lambda x: (x.free_energy))
        return cls(
            objs[0].id,
            sorted_results,
            stderr=_flatten([x.stderr for x in objs]),
            process_type=objs[0].process_type,
        )


@dataclass
class AlgorithmError:
    id: str
    error: str
