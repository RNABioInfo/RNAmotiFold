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


class GenericOutput[
    T: (ResultMFE, ResultPFC, ResultAlignment, ResultError)
](ABC):
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
    def merge(
        cls, items: list["GenericOutput[T]"]
    ) -> list["GenericOutput[T]"]:
        """Merge outputs and return a list of outputs."""
        ...

    def __repr__(self) -> str:
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(len(v)):
            together.append(f"{k[i]}={k[i]!r}")
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
            print(result_obj.tsv)
        return True

    def __str__(self) -> str:
        header = self.results[0].header
        body = "\n".join(x.tsv for x in self.results)
        return f"{header}\n{body}"

    @abstractmethod
    def to_df(self) -> DataFrame:
        """Turns Algorithm output into pandas dataframes for easy analysis and integration"""


class MFEOutput(GenericOutput[ResultMFE]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultMFE],
        stderr: list[str],
    ) -> None:
        if isinstance(result, list):
            super().__init__(
                name=name, results=sorted(result), stderr=stderr
            )
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
    def merge(
        cls: type[GenericOutput[ResultMFE]],
        items: list[GenericOutput[ResultMFE]],
    ) -> list[GenericOutput[ResultMFE]]:
        """Quick merge function for a list of algorithm outputs, no checks are built in whether they all have the same ID or anything so be careful what you input"""
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
            pandas
            .read_csv(StringIO(str(self)), sep=BaseResult.separator)
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
            super().__init__(
                name=name, results=sorted(result), stderr=stderr
            )
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
        cls: type[GenericOutput[ResultPFC]],
        items: list[GenericOutput[ResultPFC]],
    ) -> list[GenericOutput[ResultPFC]]:
        """PFC postprocessing for Single Motif mode pfc predictions filtering cases where only one motif was found"""
        returnlist: list[GenericOutput[ResultPFC]] = []
        checklist: list[str] = []
        for output in items:
            if str(output) not in checklist and len(output.results) > 1:
                checklist.append(str(output))
                returnlist.append(output)
        return returnlist

    def to_df(self) -> DataFrame:
        return (
            pandas
            .read_csv(StringIO(str(self)), sep=BaseResult.separator)
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
            super().__init__(
                name=name, results=sorted(result), stderr=stderr
            )
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
    def merge(
        cls: type[GenericOutput[ResultAlignment]],
        items: list[GenericOutput[ResultAlignment]],
    ) -> list[GenericOutput[ResultAlignment]]:
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
            pandas
            .read_csv(StringIO(str(self)), sep=BaseResult.separator)
            .sort_values(by=["id", "overall_score"])
            .reindex()
        )


class ErrorOutput(GenericOutput[ResultError]):
    def __init__(
        self,
        name: str,
        result: str | list[ResultError],
        stderr: list[str],
    ) -> None:
        if isinstance(result, list):
            super().__init__(name=name, results=result, stderr=stderr)
        else:
            relist_error: list[ResultError] = []
            relist_error.append(ResultError.from_string(name, result))
            super().__init__(name=name, results=relist_error, stderr=stderr)

    @classmethod
    def merge(
        cls: type[GenericOutput[ResultError]],
        items: list[GenericOutput[ResultError]],
    ) -> list[GenericOutput[ResultError]]:
        return list(items)

    def write_results(
        self, initiated: bool = False, destination: IO[str] = sys.stderr
    ) -> Literal[True]:
        if not initiated:
            print(
                f"Encountered issue during processing of {self.id}:",
                end="",
                file=destination,
            )
        for i in range(len(self.results)):
            print(f"{self.stderr[i]}: {self.results[i]}", file=destination)
        return True

    def to_df(self) -> DataFrame:
        return pandas.read_csv(
            StringIO(str(self)), sep=BaseResult.separator
        )


@dataclass
class AlgorithmOutput:
    results: list[MFEOutput | AlignmentOutput | PFCOutput]
    errors: list[ErrorOutput]

    @classmethod
    def from_single(
        cls,
        result: MFEOutput | ErrorOutput | AlignmentOutput | PFCOutput,
    ) -> Self:
        if isinstance(result, ErrorOutput):
            return cls(results=[], errors=[result])
        return cls(results=[result], errors=[])

    @classmethod
    def from_list(
        cls,
        result_list: (
            list[MFEOutput | ErrorOutput]
            | list[AlignmentOutput | ErrorOutput]
            | list[PFCOutput | ErrorOutput]
        ),
    ) -> Self:
        only_results: list[MFEOutput | AlignmentOutput | PFCOutput] = [
            x for x in result_list if not isinstance(x, ErrorOutput)
        ]
        only_errors: list[ErrorOutput] = [
            x for x in result_list if isinstance(x, ErrorOutput)
        ]
        return cls(results=only_results, errors=only_errors)

    def to_df(self) -> list[pandas.DataFrame] | pandas.DataFrame:
        if not self.results:
            raise ValueError(
                "Results are empty in AlgorithmOutput, cant create DataFrame"
            )
        output = type(self.results[0])
        if all(type(x) is output for x in self.results):
            merged = output.merge(self.results)  # type: ignore This throws a type error due to the generics not lining upwith the Output types, it works and is probably a lot more work that it is worth to fix, especially since its just a typing error.
            if len(merged) > 1:
                return [x.to_df().drop_duplicates() for x in merged]
            return merged[0].to_df().drop_duplicates()
        raise ValueError(
            "Different algorithms in the same AlgorithmOutput detected, can't merge"
        )
