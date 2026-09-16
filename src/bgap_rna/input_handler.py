from typing import NamedTuple, Literal, Any
from _collections_abc import Generator
from pathlib import Path
import Bio.SeqIO.FastaIO, Bio.SeqIO, Bio.SeqIO.QualityIO
import Bio.Align
from Bio.SeqRecord import SeqRecord
import gzip
import glob

#List flattening
def flatten(xss:list[list[Any]]) -> list[Any]:
    '''
    Used to flatten a lists of lists into a single list
    '''
    return [x for xs in xss for x in xs]

class algorithm_input(NamedTuple):
    id: str
    input_str: str
    call:str = ""
    
    
    def runtime_call(self) -> str:
        if self.call == "":
            raise ValueError("No call set yet for this algorithm input object")
        return " ".join([self.call, self.input_str])

    @classmethod
    def from_generators(cls,process_type:Literal["ali","mfe","pfc"],input_generator:Any) -> list['algorithm_input']:
        """Classmethod to create list of algroithm inputs (concated alignment sequences or individual mfe prediction sequences), depending on the process type
        given to the input handler by the bgap object.
        """
        returnlist:list[algorithm_input] = []
        match process_type:
            case "ali":
                for alignment in input_generator:
                    try:
                        seqs = [x.replace("-","_").replace("T","U").replace("t","u").upper() for x in alignment]
                        for seq in seqs:
                            if any(c not in "NAUCGTnaucgt+_-#" for c in set(seq)):
                                raise ValueError(
                                    "Input was neither a viable file path nor a viable RNA or DNA sequence"
                                )
                    except:
                        raise TypeError("Couldn't convert sequences from DNA to RNA")
                    concat = "#".join(seqs)
                    returnlist.append(algorithm_input("N/A",concat))
            case "mfe"|"pfc":
                for sequence in input_generator:
                    returnlist.append(algorithm_input(id=sequence.id,input_str=str(sequence.seq)))
            case _:
                raise ValueError("Unknown state detected during algorithm input generation")
        return returnlist
            
class input_handler:
    """Class to manage inputs, takes the chosen algorithm, the given inputs etc into account and creates calls from them"""
    
    def __iter__(self):
        return self

    def __next__(self):
        if self._index < len(self.user_input):
            item = self.user_input[self._index]
            self._index += 1
            return item
        else:
            self._index = 0
            raise StopIteration

    def __init__(self, process_type: Literal["mfe", "pfc", "ali"], user_input: str) -> None:
        self.user_input:list[algorithm_input] = input_handler._read_input(process_type, user_input)
        self.process_type = process_type
        self._index = 0

    @staticmethod
    def _read_input(
        process_type: Literal["mfe", "pfc", "ali"], user_input: str, id: str = "N/A"
    ) -> list[algorithm_input]:
        """Read input given by user and tries to determine if it's a file, a directory or just a single sequence.
        I really wanted to strictly type this but the number of different possible iterators make it basically impossible.
        This will always return a list of Iterators, even if there is only on MSA or seq in the file, we take care of that after."
        """
        try:
            pathd = Path(
                user_input.strip()
            )  # If this one doesn't work then it's also not a directory
            if pathd.resolve().is_file():
                result =  [input_handler._read_input_file(pathd, process_type)]
                return flatten([algorithm_input.from_generators(process_type,x) for x in result])
            elif pathd.resolve().is_dir():
                result =  [
                    input_handler._read_input_file(pathd / Path(file), process_type)
                    for file in glob.glob("*",root_dir=pathd)
                ]
                return flatten([algorithm_input.from_generators(process_type,x) for x in result])
            if any(c not in "NAUCGTnaucgt+_-#" for c in set(user_input.strip())):
                raise ValueError(
                    "Input was neither a viable file path nor a viable RNA or DNA sequence"
                )
        except OSError:
            pass
        except ValueError as e:
            raise e
          # No error was raised so we can assume that it is not a directory or file in the system but it is a viable sequence
        return [algorithm_input(id=id,input_str=user_input.strip())]
        raise ValueError("We weren't meant to get here? There was an issue with reading your input.")

    # Read input file
    @staticmethod
    def _read_input_file(
        file_path: Path, process_type: Literal["mfe", "pfc", "ali"]
    ) -> (
        Bio.SeqIO.FastaIO.FastaIterator
        | Bio.SeqIO.QualityIO.FastqPhredIterator
        | Generator[SeqRecord, None, None]
    ):
        zipped, filetype = input_handler._find_filetype(file_path)
        if (
            process_type == "ali"
        ):  # Either parse it with the Align parser if we're doing alignments or with the SeqIO if we're doing MFE/PFC calculations
            parse_func = Bio.Align.parse  # type: ignore
        else:
            parse_func = Bio.SeqIO.parse  # type: ignore
        if not zipped:
            for option in filetype:
                try:
                    return parse_func(file_path, option)  # type: ignore Both of these ignores are because of funky typing on parse from Bio
                except :
                    pass
        else:
            with gzip.open(file_path, "rt") as handle:
                for option in filetype:
                    try:
                        return parse_func(handle, option)  # type: ignore
                    except:
                        pass
        raise NotImplementedError(
            "Could not read given input file, probably because the given file type is not implemented. Try a different one or write up an issue on GitHub!"
        )

    # Finds File type based on file ending
    @staticmethod
    def _find_filetype(file_path: Path) -> tuple[bool, list[str]]:
        """Function for determining file type and zip status of user input, if file suffix is unknown this will try to use the suffix as a argument for BioSeqIO parse"""
        if file_path.suffixes[-1] == ".gz" or file_path.suffixes[-1] == ".zip":
            file_extension = file_path.suffixes[-2]
            input_zipped = True
        else:
            file_extension = file_path.suffixes[-1]
            input_zipped = False
        match file_extension:
            case (
                ".fasta"
                | ".fas"
                | ".fa"
                | ".fna"
                | ".ffn"
                | ".faa"
                | ".mpfa"
                | ".frn"
                | ".txt"
                | ".fsa"
            ):
                filetype = ["fasta", "emboss", "tabular"]
            case ".fastq" | ".fq":
                filetype = ["fastq"]
            case ".stk" | ".stockholm" | ".sto":
                filetype = ["stockholm"]
            case "phy":
                filetype = ["phylip"]
            case _:
                filetype = [file_extension]
        return (input_zipped, filetype)


if __name__ == "__main__":
    new_handler = input_handler("ali", "/home/ubuntu/ztest/test.txt")
    print(new_handler.user_input)
    print("___")
    new_handler = input_handler("ali", "/home/ubuntu/ztest")
    print(new_handler.user_input)
    print("___")
    new_handler = input_handler("mfe", "/home/ubuntu/ztest/test.txt")
    print(new_handler.user_input)
    print("___")
    new_handler = input_handler("mfe", "/home/ubuntu/ztest")
    print(new_handler.user_input)
    #new_handler = input_handler("ali", "/home/ubuntu/ztest/")
    #new_handler = input_handler("ali", "GGGGAGACCC")
    #OK we can read our input now, should work fine for now. Next step: Convert inputs into a unified format: A list of strings, each is an individual input
    #Be it an alignment, or individual sequences. Generating the right call for each input is left to the next part. This makes it a clean break between process
    #steps and is much easier to follow. Input handler reads inputs, processes them and returns a list of individual inputs to process. Since the main process knows
    #what computations it's gonna do, we can just use that. All functions here also work for the interactive mode, though it is a lot of overhead with all the processing.
    #If someone puts in a single sequence this always does the full check, though then again it is kinda worth it.
    #Actually dont return a string, for each input we return one instance of the dataclass defined above, so we can assign an id to each input and use that later during outputs!
