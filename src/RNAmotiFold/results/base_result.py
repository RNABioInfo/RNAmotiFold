import logging
from typing import Literal,Any


logger = logging.getLogger("results")

# List flattening
def flatten(xss:list[list[Any]]) -> list[Any]:
    '''
    Used to flatten a lists of lists into a single list
    '''
    return [x for xs in xss for x in xs]


class Result:
    separator: str = "\t"

    # Not yet sure how to handle motif_type it really is only interesting for single motif mode to differentiate between Internal and Bulge Loop C-Loops
    motif_type:Literal["hairpin","internal","bulge","all"] = "all"

    def __init__(self,id:str,classifier:str) -> None:
        self.id = id
        if len(classifier) == 0:
            self.classifier = "_"
        else:
            self.classifier = classifier

    def __str__(self) -> str:
        return self.tsv

    @property
    def tsv(self) -> str:
        """Returns tsv string of itself"""
        if Result.separator == r"\t":
            Result.separator = "\t"
        return Result.separator.join([str(self.__dict__[x]) for x in self.__dict__])

    @property
    def header(self) -> str:
        if Result.separator == r"\t":
            Result.separator = "\t"
        """Returns header string of itself, adapted to currently set algorithm"""
        return Result.separator.join(self.__dict__.keys())
