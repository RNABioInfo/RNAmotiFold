from RNAmotiFold.results.results import (
    ResultPFC,
    ResultMFE,
    ResultAlignment,
    ResultError,
    BaseResult,
)
from collections.abc import Iterator
from typing import Generic, TypeVar,Literal, Any,Sequence,Self
import sys
from io import StringIO
import pandas
import logging
from abc import abstractmethod,ABC
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# Typing for AlgorithmOutputs:
O = TypeVar("O",bound="GenericOutput")
T = TypeVar("T", ResultMFE, ResultPFC, ResultAlignment, ResultError)

# List flattening
def _flatten(xss: list[list[Any]]) -> list[Any]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]


class GenericOutput(ABC,Generic[T]):
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

    @classmethod
    @abstractmethod
    def merge(cls:type[O], items: list[O]) -> list[O]:
        """Merge outputs and return a list of outputs."""
        ...

    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(0, len(v)):
            together.append(
                "{key}={value!r}".format(key=k[i], value=v[i])
            )
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
    def write_results(self, initiated: bool = False) -> Literal[True]:
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
        return "\n".join([header, body])

    def to_df(self):
        return pandas.read_csv(
            StringIO(str(self)), sep=BaseResult.separator
        )


class MFEOutput(GenericOutput[ResultMFE]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultMFE],
        stderr: list[str],
    ):
        if isinstance(result, list):
            super().__init__(name=name, results=sorted(result), stderr=stderr)
        else:
            split: list[str] = result.strip().split("\n")
            reslist_mfe: list[ResultMFE] = []
            for output in split:
                res_mfe: ResultMFE = ResultMFE.from_string(
                    name, result_string=output
                )
                reslist_mfe.append(res_mfe)
            reslist_mfe.sort()
            super().__init__(
                name=name, results=sorted(reslist_mfe), stderr=stderr
            )

    @classmethod
    def merge(cls:type[O], items: list[O]) -> list[O]:
        """Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input"""
        logger.debug(f"Merging results for {items[0].id}")
        result_set: set[ResultMFE] = set()
        for obj in items:
            for res in obj.results:
                result_set.add(res)
        result_list: list[ResultMFE] = list(result_set)
        return [cls(
            items[0].id,
            result_list,
            stderr=_flatten([x.stderr for x in items]),
        )]


class PFCOutput(GenericOutput[ResultPFC]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultPFC],
        stderr: list[str],
    ):
        if isinstance(result, list):
            super().__init__(name=name, results=sorted(result), stderr=stderr)
        else:
            split = result.strip().split("\n")
            reslist_pfc: list[ResultPFC] = []
            pfc_sum = float(
                sum([float(x.split("|")[1]) for x in split])
            )
            for output in split:
                res_pfc: ResultPFC = ResultPFC.from_string(
                    id=name, result_string=output, pfc_sum=pfc_sum
                )
                reslist_pfc.append(res_pfc)
            reslist_pfc.sort(
                reverse=True
            )  # pfc has to be flipped because bigger pfc  --> more probable
            super().__init__(
                name=name, results=sorted(reslist_pfc), stderr=stderr
            )

    @classmethod
    def merge(cls:type[O], items: list[O]) -> list[O]:
        """PFC postprocessing for Single Motif mode pfc predictions filtering cases where only one motif was found"""
        returnlist: list[O] = []
        checklist: list[str] = []
        for output in items:
            if str(output) not in checklist and len(output.results) > 1:
                checklist.append(str(output))
                returnlist.append(output)
        return returnlist


class AlignmentOutput(GenericOutput[ResultAlignment]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultAlignment],
        stderr: list[str],
    ):
        if isinstance(result, list):
            super().__init__(name=name, results=sorted(result), stderr=stderr)
        else:
            split = result.strip().split("\n")
            reslist_ali: list[ResultAlignment] = []
            for output in split:
                res_ali: ResultAlignment = ResultAlignment.from_string(
                    id=name, results_string=output
                )
                reslist_ali.append(res_ali)
            reslist_ali.sort()
            super().__init__(
                name=name, results=sorted(reslist_ali), stderr=stderr
            )

    @classmethod
    def merge(cls:type[O], items: list[O]) -> list[O]:
        result_set: set[ResultAlignment] = set()
        for obj in items:
            for res in obj.results:
                result_set.add(res)
        result_list: list[ResultAlignment] = list(result_set)
        return [cls(
            items[0].id,
            result_list,
            stderr=_flatten([x.stderr for x in items]),
        )]


class ErrorOutput(GenericOutput[ResultError]):

    def __init__(
        self,
        name: str,
        result: str | list[ResultError],
        stderr: list[str],
    ):
        if isinstance(result, list):
            return super().__init__(
                name=name, results=result, stderr=stderr
            )
        else:
            relist_error: list[ResultError] = []
            relist_error.append(ResultError.from_string(name, result))
            super().__init__(
                name=name, results=relist_error, stderr=stderr
            )

    @classmethod
    def merge(cls, items: list[O]) -> list[O]:
        return list(items)


    def write_results(self, initiated: bool=False) -> Literal[True]:
        if not initiated:
            sys.stderr.write(
                f"Encountered issue during processing of {self.id}:"
            )
        for i in range(len(self.results)):
            sys.stderr.write(f"{self.stderr[i]}: {self.results[i]}")
        return True

T2 = TypeVar("T2", bound=GenericOutput)

@dataclass
class AlgorithmOutput(Generic[T2]):
    results: list[T2]
    errors: list[T2]

    @classmethod
    def from_list(
        cls,
        result_list: list[T2],
    ):
        only_results: list[T2] = [x for x in result_list if not isinstance(x, ErrorOutput)]
        only_errors: list[T2] = [x for x in result_list if isinstance(x, ErrorOutput)]
        return cls(results=only_results, errors=only_errors)

    def to_df(self):
        if not self.results:
            raise ValueError("Results are empty in AlgorithmOutput, cant create DataFrame")
        output = type(self.results[0])
        if not all(type(x) is output for x in self.results):
            raise ValueError("Different algorithms in the same AlgorithmOutput detected, can't merge")
        merged = output.merge(self.results)
        return [x.to_df().drop_duplicates() for x in merged]