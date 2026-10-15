import logging
import re
from typing import NamedTuple
from dataclasses import dataclass

logger = logging.getLogger("results")


class BaseResult:
    separator: str = "\t"

    def __init__(self, id: str, classifier: str) -> None:
        self.id = id
        if len(classifier) == 0:
            self.classifier = "_"
        else:
            self.classifier = classifier.removesuffix(
                ","
            )  # In case we use RNAmotiCes, else this does nothing

    def __str__(self) -> str:
        return self.tsv

    @property
    def tsv(self) -> str:
        """Returns tsv string of itself"""
        if BaseResult.separator == r"\t":
            BaseResult.separator = "\t"
        return BaseResult.separator.join(
            [str(self.__dict__[x]) for x in self.__dict__]
        )

    @property
    def header(self) -> str:
        if BaseResult.separator == r"\t":
            BaseResult.separator = "\t"
        """Returns header string of itself, adapted to currently set algorithm"""
        return BaseResult.separator.join(self.__dict__.keys())

    @property
    def prediction_id(self):
        return self._prediction_id

    @prediction_id.setter
    def prediction_id(self, new_id: str):
        self._prediction_id = str(new_id)


class ResultPFC(BaseResult):
    """Subclass of result for PFC results, adds pfc_value and pfc_sum attributes as well as probability calculation. Implements comparison based on pfc values."""

    def __init__(
        self,
        id: str,
        classifier: str,
        pfc_value: str | int | float,
        probability: float,
    ) -> None:
        super().__init__(id, classifier)
        self.pfc_value = float(pfc_value)
        self.probability = probability

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ResultPFC):
            return (
                self.pfc_value == other.pfc_value
                and self.classifier == other.classifier
            )
        else:
            raise NotImplementedError(
                f"Cannot compare result_pfc with {type(other)}"
            )

    def __ne__(self, other: object) -> bool:
        if isinstance(other, ResultPFC):
            return not self.__eq__(other)
        else:
            raise NotImplementedError(
                f"Cannot compare result_pfc with {type(other)}"
            )

    def __lt__(self, other: "ResultPFC") -> bool:
        return self.pfc_value < other.pfc_value

    def __le__(self, other: "ResultPFC") -> bool:
        return self.pfc_value <= other.pfc_value

    def __gt__(self, other: "ResultPFC") -> bool:
        return self.pfc_value > other.pfc_value

    def __ge__(self, other: "ResultPFC") -> bool:
        return self.pfc_value >= other.pfc_value

    @classmethod
    def from_string(
        cls, id: str, result_string: str, pfc_sum: float
    ) -> "ResultPFC":
        split_result = result_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(
            id=id,
            classifier=split_stripped_results[0],
            pfc_value=split_stripped_results[1],
            probability=round(
                float(split_stripped_results[1]) / pfc_sum, 4
            ),
        )


class ResultMFE(BaseResult):
    """Subclass of result for mfe results, has extra attributes for free energy and motBracket structure. Also implements comparison and hashing for mering structures in single motif mode"""

    def __init__(
        self,
        id: str,
        classifier: str,
        free_energy: str,
        mot_bracket: str,
    ) -> None:
        super().__init__(id, classifier)
        self.free_energy: float = float(int(free_energy) / 100)
        self.motBracket: str = mot_bracket

    # Special Dunder Method for hashing and comparing mfe results, used for merging structure in single motif mode
    def __eq__(self, other: object) -> bool:
        if isinstance(other, ResultMFE):
            return (
                self.free_energy == other.free_energy
                and self.motBracket == other.motBracket
                and self.classifier == other.classifier
            )
        else:
            raise NotImplementedError(
                f"Cannot compare result_mfe with {type(other)}"
            )

    def __ne__(self, other: object) -> bool:
        if isinstance(other, ResultMFE):
            return not self.__eq__(other)
        else:
            raise NotImplementedError(
                f"Cannot compare result_mfe with {type(other)}"
            )

    def __lt__(self, other: "ResultMFE") -> bool:
        return self.free_energy < other.free_energy

    def __le__(self, other: "ResultMFE") -> bool:
        return self.free_energy <= other.free_energy

    def __gt__(self, other: "ResultMFE") -> bool:
        return self.free_energy > other.free_energy

    def __ge__(self, other: "ResultMFE") -> bool:
        return self.free_energy >= other.free_energy

    def __hash__(self) -> int:
        return hash(
            (self.free_energy, self.motBracket, self.classifier)
        )

    @property
    def dot_bracket(self) -> str:
        return self._dot_bracket

    @dot_bracket.setter
    def dot_bracket(self, structure_string: str) -> None:
        for c in set(self.classifier):
            structure_string = structure_string.replace(c, ".")
        self._dot_bracket = structure_string

    @classmethod
    def from_string(cls, id: str, result_string: str) -> "ResultMFE":
        split_result = result_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(
            id=id,
            classifier=split_stripped_results[0],
            free_energy=split_stripped_results[1],
            mot_bracket=split_stripped_results[2],
        )


class ResultAlignment(BaseResult):
    """Subclass of result for alignment results, adds alignment_score and motBracket attributes as well as special implementation for tsv and header functions"""

    def __init__(
        self, id: str, classifier: str, score: str, motBracket: str
    ) -> None:
        super().__init__(id, classifier)
        split_score: alignment_score = alignment_score.from_string(
            score
        )
        self.energy: float = split_score.energy
        self.covariance: float = split_score.covariance
        self.motif_score: float = split_score.motif
        self.overall_score: float = split_score.overall
        self.motBracket: str = motBracket

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ResultAlignment):
            return (
                self.energy == other.energy
                and self.covariance == other.covariance
                and self.overall_score == other.overall_score
                and self.motif_score == other.motif_score
                and self.motBracket == other.motBracket
                and self.classifier == other.classifier
            )
        else:
            raise NotImplementedError(
                f"Cannot compare result_alignment with {type(other)}"
            )

    def __ne__(self, other: object) -> bool:
        if isinstance(other, ResultAlignment):
            return not self.__eq__(other)
        else:
            raise NotImplementedError(
                f"Cannot compare result_alignment with {type(other)}"
            )

    def __lt__(self, other: "ResultAlignment") -> bool:
        return self.overall_score < other.overall_score

    def __le__(self, other: "ResultAlignment") -> bool:
        return self.overall_score <= other.overall_score

    def __gt__(self, other: "ResultAlignment") -> bool:
        return self.overall_score > other.overall_score

    def __ge__(self, other: "ResultAlignment") -> bool:
        return self.overall_score >= other.overall_score

    @classmethod
    def from_string(
        cls, id: str, results_string: str
    ) -> "ResultAlignment":
        split_result = results_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(
            id=id,
            classifier=split_stripped_results[0],
            score=split_stripped_results[1],
            motBracket=split_stripped_results[2],
        )

    def __hash__(self) -> int:
        return hash((
            self.energy,
            self.covariance,
            self.overall_score,
            self.motif_score,
            self.motBracket,
            self.classifier,
        ))



class alignment_score(NamedTuple):
    energy: float
    covariance: float
    motif: float
    overall: float

    @classmethod
    def from_string(cls, score_string: str):
        nums = [
            float(x)
            for x in re.findall(r"-?\d+(?:\.\d+)?", score_string)
        ]
        return cls(
            overall=nums[0],
            energy=nums[1],
            covariance=nums[2],
            motif=nums[3],
        )

    def to_string(self, sep: str):
        return sep.join(
            [
                str(self.overall),
                str(self.energy),
                str(self.covariance),
                str(self.motif),
            ]
        )


@dataclass
class ResultError(BaseResult):
    def __init__(self, id: str, errstr: str):
        super().__init__(id=id, classifier="ERROR")
        self.error = errstr

    @classmethod
    def from_string(cls, id: str, errstr: str):
        return cls(id=id, errstr=errstr)
