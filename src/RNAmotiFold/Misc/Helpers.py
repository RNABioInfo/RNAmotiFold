from numpy import float128


# List flattening
def flatten[T](xss: list[list[T]]) -> list[T]:
    """
    Used to flatten a lists of lists into a single list
    """

    return [x for xs in xss for x in xs]

def interval_check[T: (float, float128)](
    LowerBound: T, UpperBound: T, Value: T
) -> bool:
    return bool(LowerBound <= Value <= UpperBound)