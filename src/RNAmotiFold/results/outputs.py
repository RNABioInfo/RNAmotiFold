from RNAmotiFold.results.results import ResultPFC, ResultMFE, ResultAlignment, ResultError, BaseResult
from collections.abc import Iterator
from typing import Generic, TypeVar
from typing import Literal, Any
import sys
from io import StringIO
import pandas
import logging
from typing import Self


logger = logging.getLogger(__name__)

# Typing for AlgorithmOutputs:
T = TypeVar("T", ResultMFE, ResultPFC, ResultAlignment, ResultError)


class AlgorithmOutput(Generic[T]):
    """Generic AlgorithmOutput class to inherit from with MFE/PFC/Alignment Subclasses"""

    def __iter__(self) -> Iterator[T]:
        return self

    def __next__(self) -> T:
        if self._index < len(self.results):
            item = self.results[self._index]
            self._index += 1
            return item

        self._index = 0
        raise StopIteration

    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(0, len(v)):
            together.append("{key}={value!r}".format(key=k[i], value=v[i]))
        return f"{classname}({', '.join(together)})"

    def __init__(
        self,
        name: str,
        results: list[T],
        stderr: list[str],
        motif: Literal["hairpin", "internal", "bulge", "all"] = "all",
    ) -> None:
        self.id: str = name
        self.results: list[T] = results
        self.stderr = stderr
        self.motif_type = motif
        self._index = 0

    @property
    def prediction_id(self) -> str:
        return self._prediction_id
    
    @prediction_id.setter
    def prediction_id(self,new_id:int):
        for result in self.results:
            result.prediction_id = str(new_id)
        self._prediction_id = str(new_id)

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

    def __str__(self) -> str:
        header = self.results[0].header
        body = "\n".join(x.tsv for x in self.results)
        return "\n".join([header,body])

    def to_df(self):
        return pandas.read_csv(StringIO(str(self)),sep=BaseResult.separator)

class MFEOutput(AlgorithmOutput[ResultMFE]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultMFE],
        stderr: list[str],
    ):
        if isinstance(result, list):
            super().__init__(name=name, results=result, stderr=stderr)
        else:
            split: list[str] = result.strip().split("\n")
            reslist_mfe: list[ResultMFE] = []
            for output in split:
                res_mfe: ResultMFE = ResultMFE.from_string(name, result_string=output)
                reslist_mfe.append(res_mfe)
            reslist_mfe.sort()
            super().__init__(name=name, results=reslist_mfe, stderr=stderr)

    @classmethod
    def merge_mfe(cls, objs: list[Self]) -> Self:
        """Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input"""
        logger.debug(f"Merging results for {objs[0].id}")
        result_set: set[ResultMFE] = set()
        for obj in objs:
            for res in obj.results:
                result_set.add(res)
        result_list: list[ResultMFE] = list(result_set)
        return cls(objs[0].id, result_list, stderr=_flatten([x.stderr for x in objs]))

class PFCOutput(AlgorithmOutput[ResultPFC]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultPFC],
        stderr: list[str],
    ):
        if isinstance(result, list):
            super().__init__(name=name, results=result, stderr=stderr)
        else:
            split = result.strip().split("\n")
            reslist_pfc: list[ResultPFC] = []
            pfc_sum = float(sum([float(x.split("|")[1]) for x in split]))
            for output in split:
                res_pfc: ResultPFC = ResultPFC.from_string(
                    id=name, result_string=output, pfc_sum=pfc_sum
                )
                reslist_pfc.append(res_pfc)
            reslist_pfc.sort(
                reverse=True
            )  # pfc has to be flipped because bigger pfc  --> more probable
            super().__init__(name=name, results=reslist_pfc, stderr=stderr)

    @staticmethod
    def merge_pfc(objs: list["PFCOutput"]) -> list["PFCOutput"]:
        """PFC postprocessing for Single Motif mode pfc predictions filtering cases where only one motif was found"""
        returnlist: list[PFCOutput] = []
        checklist: list[str] = []
        for output in objs:
            if str(output) not in checklist and len(output.results) > 1:
                checklist.append(str(output))
                returnlist.append(output)
        return returnlist


class AlignmentOutput(AlgorithmOutput[ResultAlignment]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultAlignment],
        stderr: list[str],
    ):
        if isinstance(result, list):
            super().__init__(name=name, results=result, stderr=stderr)
        else:
            split = result.strip().split("\n")
            reslist_ali: list[ResultAlignment] = []
            for output in split:
                res_ali: ResultAlignment = ResultAlignment.from_string(
                    id=name, results_string=output
                )
                reslist_ali.append(res_ali)
            reslist_ali.sort()
            super().__init__(name=name, results=reslist_ali, stderr=stderr)

    @classmethod
    def merge_ali(cls, objs: list[Self]) -> Self:
        result_set: set[ResultAlignment] = set()
        for obj in objs:
            for res in obj.results:
                result_set.add(res)
        result_list: list[ResultAlignment] = list(result_set)
        return cls(objs[0].id, result_list, stderr=_flatten([x.stderr for x in objs]))



# List flattening
def _flatten(xss: list[list[Any]]) -> list[Any]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]


class ErrorOutput(AlgorithmOutput[ResultError]):

    def __init__(self, name: str, result: str | list[ResultError], stderr: list[str]):
        if isinstance(result, list):
            return super().__init__(name=name, results=result, stderr=stderr)
        else:
            relist_error: list[ResultError] = []
            relist_error.append(ResultError.from_string(name, result))
            super().__init__(name=name, results=relist_error, stderr=stderr)

    def write_results(self, initiated: bool) -> Literal[True]:
        if not initiated:
            sys.stderr.write(f"Encountered issue during processing of {self.id}:")
        for i in range(len(self.results)):
            sys.stderr.write(f"{self.stderr[i]}: {self.results[i]}")
        return True
