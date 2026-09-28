from RNAmotiFold.results.base_result import Result

class ResultPFC(Result):
    """Subclass of result for PFC results, adds pfc_value and pfc_sum attributes as well as probability calculation. Implements comparison based on pfc values."""

    def __init__(
        self, id: str, classifier: str, pfc_value: str | int | float, probability: float
    ) -> None:
        super().__init__(id, classifier)
        self.pfc_value = float(pfc_value)
        self.probability = probability

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ResultPFC):
            return self.pfc_value == other.pfc_value and self.classifier == other.classifier
        else:
            raise NotImplementedError(f"Cannot compare result_pfc with {type(other)}")

    def __ne__(self, other: object) -> bool:
        if isinstance(other, ResultPFC):
            return not self.__eq__(other)
        else:
            raise NotImplementedError(f"Cannot compare result_pfc with {type(other)}")

    def __lt__(self, other: "ResultPFC") -> bool:
        return self.pfc_value < other.pfc_value

    def __le__(self, other: "ResultPFC") -> bool:
        return self.pfc_value <= other.pfc_value

    def __gt__(self, other: "ResultPFC") -> bool:
        return self.pfc_value > other.pfc_value

    def __ge__(self, other: "ResultPFC") -> bool:
        return self.pfc_value >= other.pfc_value

    @classmethod
    def from_string(cls, id: str, result_string: str, pfc_sum: float) -> "ResultPFC":
        split_result = result_string.strip().split("|")
        split_stripped_results = [x.strip() for x in split_result]
        return cls(
            id=id,
            classifier=split_stripped_results[0],
            pfc_value=split_stripped_results[1],
            probability=round(float(split_stripped_results[1]) / pfc_sum, 4),
        )
