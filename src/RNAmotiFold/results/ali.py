import re
from typing import NamedTuple
import src.RNAmotiFold.results.base_result


class alignment_score(NamedTuple):
    energy: float
    covariance: float
    motif: float
    overall: float

    @classmethod
    def from_string(cls, score_string: str):
        nums = [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?", score_string)]
        return cls(overall=nums[0], energy=nums[1], covariance=nums[2], motif=nums[3])

    def to_string(self, sep: str):
        return sep.join(
            [str(self.overall), str(self.energy), str(self.covariance), str(self.motif)]
        )


class result_alignment(src.RNAmotiFold.results.base_result.result):
    """Subclass of result for alignment results, adds alignment_score and motBracket attributes as well as special implementation for tsv and header functions"""

    def __init__(self, id: str, classifier: str, score: str, motBracket: str) -> None:
        super().__init__(id, classifier)
        self.score: alignment_score = alignment_score.from_string(score)
        self.motBracket: str = motBracket

    def __eq__(self, other: object) -> bool:
        if isinstance(other, result_alignment):
            return (
                self.score == other.score
                and self.motBracket == other.motBracket
                and self.classifier == other.classifier
            )
        else:
            raise NotImplementedError(f"Cannot compare result_alignment with {type(other)}")

    def __ne__(self, other: object) -> bool:
        if isinstance(other, result_alignment):
            return not self.__eq__(other)
        else:
            raise NotImplementedError(f"Cannot compare result_alignment with {type(other)}")

    def __lt__(self, other: "result_alignment") -> bool:
        return self.score.overall < other.score.overall

    def __le__(self, other: "result_alignment") -> bool:
        return self.score.overall <= other.score.overall

    def __gt__(self, other: "result_alignment") -> bool:
        return self.score.overall > other.score.overall

    def __ge__(self, other: "result_alignment") -> bool:
        return self.score.overall >= other.score.overall

    @property
    def overall_score(self) -> float:
        return self.score.overall

    @property
    def energy(self) -> float:
        return self.score.energy

    @property
    def covariance(self) -> float:
        return self.score.covariance

    @property
    def motif_score(self) -> float:
        return self.score.motif

    @property
    def header(self) -> str:
        return src.RNAmotiFold.results.base_result.result.separator.join(
            [
                "ID",
                "Motif",
                "Total Score",
                "Free Energy",
                "Covariance Score",
                "Motif Score",
                "MotBracket",
            ]
        )

    @property
    def tsv(self):
        return src.RNAmotiFold.results.base_result.result.separator.join(
            [
                self.id,
                self.classifier,
                str(self.score.overall),
                str(self.score.energy),
                str(self.score.covariance),
                str(self.score.motif),
                self.motBracket,
            ]
        )

    @classmethod
    def from_string(cls, id: str, results_string: str) -> "result_alignment":
        split_result = results_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(
            id=id,
            classifier=split_stripped_results[0],
            score=split_stripped_results[1],
            motBracket=split_stripped_results[2],
        )
