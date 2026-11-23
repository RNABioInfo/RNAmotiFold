from typing import TypeVar

from numpy import float128

T = TypeVar("T")
IntervalValue = TypeVar("IntervalValue", float, float128)


# List flattening
def flatten(xss: list[list[T]]) -> list[T]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]


def interval_check(
    LowerBound: IntervalValue, UpperBound: IntervalValue, Value: IntervalValue
) -> bool:
    return bool(LowerBound <= Value <= UpperBound)
