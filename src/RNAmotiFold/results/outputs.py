import logging
import sys
from abc import ABC, abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass
from io import StringIO
from typing import IO, Any, Literal, Self

import pandas
from pandas import DataFrame
from RNAmotiFold.results.results import (
    BaseResult,
    ResultAlignment,
    ResultError,
    ResultMFE,
    ResultPFC,
)

logger = logging.getLogger(__name__)

# Typing for AlgorithmOutputs:


# List flattening
def _flatten(xss: list[list[Any]]) -> list[Any]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]

class GenericOutput[T: (ResultMFE, ResultPFC, ResultAlignment, ResultError)](ABC):
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

    def __repr__(self) -> str:
        classname = type(self).__name__
        together = [f"{key}={value!r}" for key, value in self.__dict__.items()]
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
    def write_results(
        self, initiated: bool = False, destination: IO[str] = sys.stdout
    ) -> Literal[True]:
        """Header and results written with this function will be in csv format using the classwide results.separator variable"""
        for err in self.stderr:
            if len(err.strip()) > 0:
                logger.warning(self.id + ": " + err.strip())
        if not initiated:
            logger.debug("Starting result writing")
            print(self.results[0].header, file=destination)
        for result_obj in self.results:
            print(result_obj.tsv, file=destination)
        return True

    def __str__(self) -> str:
        header = self.results[0].header
        body = "\n".join(x.tsv for x in self.results)
        return f"{header}\n{body}"

    @abstractmethod
    def to_df(self) -> DataFrame:
        """Turns Algorithm output into pandas dataframes for easy analysis and integration"""
        ...

    @classmethod
    @abstractmethod
    def merge[subT:(MFEOutput,PFCOutput,AlignmentOutput)](cls:type[subT],items:list[subT]) -> list[subT]:
        ...


class MFEOutput(GenericOutput[ResultMFE]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultMFE],
        stderr: list[str],
    ) -> None:
        if isinstance(result, list):
            super().__init__(name=name, results=sorted(result), stderr=stderr)
        else:
            split: list[str] = result.strip().split("\n")
            reslist_mfe: list[ResultMFE] = []
            for output in split:
                res_mfe: ResultMFE = ResultMFE.from_string(name, result_string=output)
                reslist_mfe.append(res_mfe)
            reslist_mfe.sort()
            super().__init__(name=name, results=sorted(reslist_mfe), stderr=stderr)

    @classmethod
    def merge(
        cls: type["MFEOutput"],
        items: list["MFEOutput"],
    ) -> list["MFEOutput"]:
        """Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input"""
        if len(items) == 0:
            return []
        logger.debug(f"Merging results for {items[0].id}")
        result_set: set[ResultMFE] = set()
        for obj in items:
            result_set.update(obj.results)
        result_list: list[ResultMFE] = list(result_set)
        return [
            cls(
                items[0].id,
                result_list,
                stderr=_flatten([x.stderr for x in items]),
            )
        ]

    def to_df(self) -> DataFrame:
        return (
            pandas.read_csv(StringIO(str(self)), sep=BaseResult.separator)
            .sort_values(by=["id", "free_energy"])
            .reindex()
        )


class PFCOutput(GenericOutput[ResultPFC]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultPFC],
        stderr: list[str],
    ) -> None:
        if isinstance(result, list):
            super().__init__(name=name, results=sorted(result), stderr=stderr)
        else:
            split = result.strip().split("\n")
            reslist_pfc: list[ResultPFC] = []
            pfc_sum = float(sum(float(x.split("|")[1]) for x in split))
            for output in split:
                res_pfc: ResultPFC = ResultPFC.from_string(
                    id=name, result_string=output, pfc_sum=pfc_sum
                )
                reslist_pfc.append(res_pfc)
            reslist_pfc.sort(
                reverse=True
            )  # pfc has to be flipped because bigger pfc  --> more probable
            super().__init__(
                name=name,
                results=sorted(reslist_pfc, reverse=True),
                stderr=stderr,
            )

    @classmethod
    def merge(
        cls: type["PFCOutput"],
        items: list["PFCOutput"],
    ) -> list["PFCOutput"]:
        """PFC postprocessing for Single Motif mode pfc predictions filtering cases where only one motif was found"""
        if len(items) == 0:
            return []
        returnlist: list[PFCOutput] = []
        checklist: list[str] = []
        for output in items:
            if str(output) not in checklist and len(output.results) > 1:
                checklist.append(str(output))
                returnlist.append(output)
        return returnlist

    def to_df(self) -> DataFrame:
        return (
            pandas.read_csv(StringIO(str(self)), sep=BaseResult.separator)
            .sort_values(by=["id", "probability"], ascending=[True, False])
            .reindex()
        )


class AlignmentOutput(GenericOutput[ResultAlignment]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultAlignment],
        stderr: list[str],
    ) -> None:
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
            super().__init__(name=name, results=sorted(reslist_ali), stderr=stderr)

    @classmethod
    def merge(
        cls: type["AlignmentOutput"],
        items: list["AlignmentOutput"],
    ) -> list["AlignmentOutput"]:
        if len(items) == 0:
            return []
        result_set: set[ResultAlignment] = set()
        for obj in items:
            result_set.update(obj.results)
        result_list: list[ResultAlignment] = list(result_set)
        return [
            cls(
                items[0].id,
                result_list,
                stderr=_flatten([x.stderr for x in items]),
            )
        ]

    def to_df(self) -> DataFrame:
        return (
            pandas.read_csv(StringIO(str(self)), sep=BaseResult.separator)
            .sort_values(by=["id", "overall_score"])
            .reindex()
        )


class ErrorOutput(GenericOutput[ResultError]):
    def __init__(
        self,
        name: str,
        call: str,
        result: str | list[ResultError],
        stderr: list[str],
    ) -> None:
        if isinstance(result, list):
            super().__init__(name=name, results=result, stderr=stderr)
            self.call = call
        else:
            relist_error: list[ResultError] = []
            relist_error.append(ResultError.from_string(name, result))
            super().__init__(name=name, results=relist_error, stderr=stderr)
            self.call = call

    @classmethod
    def merge(
        cls: type["ErrorOutput"],
        items: list["ErrorOutput"],
    ) -> list["ErrorOutput"]:
        if len(items) == 0:
            return []
        return list(items)

    def write_results(
        self, initiated: bool = False, destination: IO[str] = sys.stderr
    ) -> Literal[True]:
        if not initiated:
            print(
                f"Encountered issue while processing {self.id} with call {self.call}:",
                end="",
                file=destination,
            )
        for i in range(len(self.results)):
            print(f"{self.stderr[i]}: {self.results[i]}", file=destination)
        return True

    def to_df(self) -> DataFrame:
        return pandas.read_csv(StringIO(str(self)), sep=BaseResult.separator)

@dataclass
class AlgorithmOutput[T: (MFEOutput,PFCOutput,AlignmentOutput)]:
    results:list[T] 
    errors :list[ErrorOutput]
    result_type:type[T]

    @classmethod
    def from_list(
        cls,
        result_list: list[T | ErrorOutput],
        result_type:type[T]
    ) -> Self:
        only_errors = [x for x in result_list if isinstance(x, ErrorOutput)]
        only_results = [x for x in result_list if isinstance(x, result_type)]
        return cls(only_results, only_errors,result_type)

    def to_df(self) -> list[pandas.DataFrame] | pandas.DataFrame:
        if len(self.results) == 0:
            raise ValueError(
                "Results are empty in AlgorithmOutput, cant create DataFrame"
            )
        merged = self.result_type.merge(self.results) #type:ignore Error caused by union can be fixed with protocol but not worth it
        return [x.to_df().drop_duplicates() for x in merged]
