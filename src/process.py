import sys
import gzip
import argparse
import subprocess
import configparser
import multiprocessing
from pathlib import Path
from typing import Iterator
from typing import Optional
from typing import Generator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from requests import get
import src.results as results
from contextlib import redirect_stdout

# A Python class for making Bellman's GAP more convenient to use
# Just create a class instances, feed it with the call arguments you need
# and it'll create a call from the arguments classified as runtime arguments.
# These are: motif_source, motif_orientation, kvalue, hishape_mode, shape_level, energy


class bgap_rna:

    def __repr__(self) -> str:
        return f"bgap_rna call {self.call}"

    def __str__(self) -> str:
        return f"{self.algorithm}:{self.__dict__}"

    @staticmethod
    def _get_cmd_arguments() -> argparse.Namespace:
        parser = argparse.ArgumentParser(
            prog="bgap_rna class",
            description="Argument parser for the bgap_rna class, allowing construction through runtime arguments of your script. Arguments have no defaults to keep in line with the idea that only set parameters are considered.",
            epilog="The helpline is currently at capacity, please hold...",
        )
        parser.add_argument(
            "-i",
            "--input",
            help="Set input path or input sequence. File formats fasta, fastq and stockholm are supported. File compression .zip and .gzip are also supported.",
            type=str,
            dest="input",
            nargs="?",
            required=True,
        )
        # Command line arguments that control which algorithm is called with which options.
        # If you add your own partition function algorithm and want the output to have probabilities be sure to add pfc at the end of the name! This tag is used to recognize partition function algorithms by the script.
        parser.add_argument(
            "-a",
            "--algorithm",
            help="Specify which algorithm should be used, prebuild choices are: motmfepretty, motpfc, motshapeX, motshapeX_pfc, mothishape, mothishape_h_pfc, mothishape_b_pfc, mothishape_m_pfc. If you want to run a mothishape you need to specify the mode with -p, if you run mothishape with pfc enter the full name (e.g. mothishape_h_pfc). Paritition function instances should be marked with pfc at the end.",
            type=str,
            nargs="?",
            dest="algorithm",
            required=True,
        )
        parser.add_argument(
            "-pa",
            "--path",
            help="Specify folder where your algorithm is compiled. Default is parent of the folder where this file is located.",
            type=Path,
            default=Path(__file__).resolve().parents[1],
            dest="algorithm_folder",
            required=True,
        )
        parser.add_argument(
            "-s",
            "--subopt",
            help="Specify if subopt folding should be used. Not useable with partition function implementations. Default is off",
            action="store_true",
            dest="subopt",
        )
        parser.add_argument(
            "-Q",
            "--motif_source",
            help="Specify from which database motifs should be used, 1 = BGSU, 2 = Rfam, 3 = both. Default is 3",
            choices=[1, 2, 3],
            type=int,
            dest="motif_source",
        )
        parser.add_argument(
            "-b",
            "--orientation",
            help="Specify motif orientation: 1 = 5'-> 3',  2 = 3' -> 5' or 3 = both. Default is 3",
            choices=[1, 2, 3],
            type=int,
            dest="motif_orientation",
        )
        parser.add_argument(
            "-k",
            "--kvalue",
            help="Specify k for k-best classes get classified. Default is k = 5",
            type=int,
            dest="kvalue",
        )
        parser.add_argument(
            "-q",
            "--shape_level",
            help="Set shape abstraction level. Default is 3",
            choices=[1, 2, 3, 4, 5],
            type=int,
            dest="shape_level",
        )
        parser.add_argument(
            "-e",
            "--energy",
            help="Set energy range in kcal/mol. Default is 10.0",
            type=float,
            dest="energy",
        )
        parser.add_argument(
            "-p",
            "--hishape",
            help="Set hishape mode, default is h",
            choices=["h", "m", "b"],
            type=str,
            dest="hishape_mode",
        )
        parser.add_argument(
            "-t",
            "--temperature",
            help="Scale energy parameters for folding to given temperature in Celsius. Default is 37",
            type=float | int,
            dest="temperature",
        )
        parser.add_argument(
            "-c",
            help="Set energy range in %. Is overwritten if -e is set. Default is 10.0",
            type=float | int,
            dest="energy_percent",
        )
        parser.add_argument(
            "-u",
            help="Allow lonely base pairs, 1 = yes, 0 = no. Default is 0",
            type=int,
            dest="basepairs",
        )
        parser.add_argument(
            "--time",
            help="Activate time logging, activating this will run predictions with unix time utility. Default is off",
            action="store_true",
            dest="time",
        )
        args = parser.parse_args()
        return args

    @staticmethod
    def _worker_funct(
        call: str,
        iq: multiprocessing.Queue,
        oq: multiprocessing.Queue,
    ) -> list:
        while True:
            record = iq.get()  # type:Optional[SeqIO.SeqRecord]
            if record is None:
                break
            else:
                if "T" in str(record.seq):
                    record.seq = record.seq.transcribe()
                subprocess_output = subprocess.run(
                    call + str(record.seq), text=True, capture_output=True, shell=True
                )
            if not subprocess_output.returncode:
                result = results.algorithm_output(
                    record.id, subprocess_output.stdout, subprocess_output.stderr
                )
            else:
                result = results.error(record.id, subprocess_output.stderr)
            oq.put(result)

    @staticmethod
    def get_current_motif_version(attempts=5) -> str:
        i = 0
        while i < attempts:
            response = get("http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json")
            if response.status_code == 200:
                return response.headers["Content-disposition"].split("=")[1].split("_")[1][:-5]
        else:
            i += 1
        raise ConnectionError(
            f"Could not establish connection to API server in {attempts} attempts."
        )

    # Compare local RNA 3D Motif Atlas version against current, input should be installed motif version findable under RNALoops/src/data/config
    @staticmethod
    def check_motif_versions(installed: str) -> bool:
        try:
            current = bgap_rna.get_current_motif_version()
        except ConnectionError as error:
            raise ConnectionError(
                "Could not connect to RNA 3D Motif Atlas server to check for current motif versions."
            )
        if current == installed:
            return False
        else:
            return True

    @classmethod
    def from_argparse(cls, cmd_args: Optional[argparse.Namespace] = None):
        """argparse alternative constructor."""
        if cmd_args is None:
            cmd_args = bgap_rna._get_cmd_arguments()
        return cls(
            cmd_args.input,
            cmd_args.algorithm,
            cmd_args.algorithm_folder,
            cmd_args.motif_source,
            cmd_args.motif_orientation,
            cmd_args.kvalue,
            cmd_args.shape_level,
            cmd_args.energy,
            cmd_args.hishape_mode,
            cmd_args.temperature,
            cmd_args.energy_percent,
            cmd_args.basepairs,
            cmd_args.subopt,
            cmd_args.time,
        )

    @classmethod
    def from_config(cls, config_path: str | Path | configparser.ConfigParser):
        """ConfigParse alternative constructor. Accepts a string/path to a config file or a ConfigParser. Unused parameters should be left empty or be set to None."""
        if isinstance(config_path, str) or isinstance(config_path, Path):
            config = configparser.ConfigParser(allow_no_value=True)
            config.read_file(open(config_path))
            return cls(
                config["PROCESS"]["input"],
                config["PROCESS"]["algorithm"],
                config["PROCESS"]["algorithm_folder"],
                config.getint("PROCESS", "motif_source"),
                config.getint("PROCESS", "motif_orientation"),
                config.getint("PROCESS", "kvalue"),
                config.getint("PROCESS", "shape_level"),
                config.getfloat("PROCESS", "energy"),
                config["PROCESS"]["hishape_mode"],
                config.getfloat("PROCESS", "temperature"),
                config.getfloat("PROCESS", "energy_percent"),
                config.getint("PROCESS", "basepairs"),
                config.getboolean("PROCESS", "subopt"),
                config.getboolean("PROCESS", "time"),
            )
        elif isinstance(config_path, configparser.ConfigParser):
            config_path.allow_no_value = True
            return cls(
                config_path["PROCESS"]["input"],
                config_path["PROCESS"]["algorithm"],
                config_path["PROCESS"]["algorithm_folder"],
                config_path.getint("PROCESS", "motif_source"),
                config_path.getint("PROCESS", "motif_orientation"),
                config_path.getint("PROCESS", "kvalue"),
                config_path.getint("PROCESS", "shape_level"),
                config_path.getfloat("PROCESS", "energy"),
                config_path["PROCESS"]["hishape_mode"],
                config_path.getfloat("PROCESS", "temperature"),
                config_path.getfloat("PROCESS", "energy_percent"),
                config_path.getint("PROCESS", "basepairs"),
                config_path.getboolean("PROCESS", "subopt"),
                config_path.getboolean("PROCESS", "time"),
            )

    def __init__(
        self,
        input: str | SeqRecord | list[str | SeqRecord],
        algorithm: str,
        algorithm_folder: str | Path,
        motif_source: Optional[int] = None,
        motif_orientation: Optional[int] = None,
        kvalue: Optional[int] = None,
        shape_level: Optional[int] = None,
        energy: Optional[float] = None,
        hishape_mode: Optional[str] = None,
        temperature: Optional[float | int] = None,
        energy_percent: Optional[float | int] = None,
        allowLonelyBasepairs: Optional[int] = None,
        subopt: bool = False,
        time: bool = False,
        pfc_filtering: bool = False,
    ):
        self.subopt = subopt
        self.time = time
        self.motif_source = motif_source
        self.motif_orientation = motif_orientation
        self.kvalue = kvalue
        self.shape_level = shape_level
        self.energy = energy
        self.hishape_mode = hishape_mode
        self.temperature = temperature
        self.algorithm = algorithm
        self.energy_percent = energy_percent
        self.allowLonelyBasepairs = allowLonelyBasepairs
        self.input = (
            input
        )  # type: SeqIO.FastaIO.FastaIterator | SeqIO.QualityIO.FastqPhredIterator | Generator[SeqRecord, None, None] | Iterator[list]
        self.algorithm_path = algorithm_folder  # Checks for the algorithm in the given folder, if it isn't there attempts to compile the algorithm
        self.pfc_filtering = pfc_filtering  # Set this to enable/disable filtering of partition function outputs to only return results with a probability over 0.0001

    # Slightly controversial addition, if custom_call is set it permanently overwrites the default call and is even returned whenever
    # the standard self.call is asked for. This avoids duplicating and overcomplicating code down the line. Deleting this will return normal calls
    @property
    def custom_call(self):
        return self._custom_call

    @custom_call.setter
    def custom_call(self, call: str):
        self._custom_call = call

    @custom_call.deleter
    def custom_call(self):
        del self._custom_call

    @property
    def temperature(self):
        return self._temperature

    @temperature.setter
    def temperature(self, temp: Optional[float]):
        if temp is None or float or int:
            self._temperature = temp
        else:
            raise ValueError("Temperature can only be None or a Float Value.")

    @property
    def energy_percent(self):
        return self._energy_perc

    @energy_percent.setter
    def energy_percent(self, range):
        if range is None or float or int:
            self._energy_perc = range
        else:
            raise ValueError("Energy range can only be a float or int.")

    @property
    def allowLonelyBasepairs(self):
        return self._allowLonelyBasepairs

    @allowLonelyBasepairs.setter
    def allowLonelyBasepairs(self, val: Optional[int]):
        if val is None or val in [0, 1]:
            self._allowLonelyBasepairs = val
        else:
            raise ValueError("allowLonelyBasepairs can only be 0 or 1.")

    # Choose motif source, 1 = RNA 3D Motif Atlas, 2 = RMFam, 3 = Both
    @property
    def motif_source(self):
        return self._motif_source

    @motif_source.setter
    def motif_source(self, Q: Optional[int]):
        if Q is None or Q in [1, 2, 3]:
            self._motif_source = Q
        else:
            raise ValueError(
                "Motif source can only be 1 = RNA 3D Motif Atlas , 2 = RMfam , 3 = Both."
            )

    # Choose motif orientation, 1 = 5'->3' only, 2 = 3'->5' only, 3= both
    @property
    def motif_orientation(self):
        return self._motif_orientation

    @motif_orientation.setter
    def motif_orientation(self, b: Optional[int]):
        if b is None or b in [1, 2, 3]:
            self._motif_orientation = b
        else:
            raise ValueError("Motif direction can only be 1 = 5'->3' , 2 = 3'->5' , 3 = Both.")

    # Set shape abstraction level, viable inputs are 1-5
    @property
    def shape_level(self):
        return self._shape_level

    @shape_level.setter
    def shape_level(self, q: Optional[int]):
        if q is None or q in [1, 2, 3, 4, 5]:
            self._shape_level = q
        else:
            raise ValueError(
                "Shape level can only be set to levels 5 (most abstract) - 1 (least abstract)."
            )

    # Set hishape abstraction level, viable inputs are h, b, m
    @property
    def hishape_mode(self):
        return self._hishape

    @hishape_mode.setter
    def hishape_mode(self, h: Optional[str]):
        if h is None or h in ["h", "m", "b"]:
            self._hishape = h
        else:
            raise ValueError("HiShape mode can only be set to h,b or m.")

    # Set energy range for suboptimal candidate calcualtions
    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, e: Optional[float | int]):
        if e is not None and e > 0:
            self._energy = float(e)
        elif e is None:
            self._energy = e
        else:
            raise ValueError("Energy range cannot be lower than 0.")

    # Set kvalue fpr kbest and kbacktracing
    @property
    def kvalue(self):
        return self._kvalue

    @kvalue.setter
    def kvalue(self, k):
        if k is not None and k > 0:
            self._kvalue = k
        elif k is None:
            self._kvalue = k
        else:
            raise ValueError("Kvalue cannot be 0 or lower.")

    # Finds path to your chosen algorithm, if it does not exist i attempts to compile the algorithm
    @property
    def algorithm_path(self):
        return self._algorithm_path

    @algorithm_path.setter
    def algorithm_path(self, home_folder: str | Path):
        self._algorithm_path = str(Path(home_folder).joinpath(self.algorithm))

    # Builds the algorithm name string, checks for _pfc at the end and builds hishape name and adds _subopt
    @property
    def algorithm(self):
        return self._algorithm

    @algorithm.setter
    def algorithm(self, alg):
        if alg == "mothishape":
            try:
                alg = alg + "_" + self.hishape_mode
            except TypeError:
                raise TypeError("Specify a hishape mode with -q if you want to use hishapes.")
        if self.subopt:
            alg = alg + "_subopt"
        if alg[-3:] == "pfc":
            self.pfc = True
            if self.subopt:
                raise ValueError("Partition function can't be used with subopt")
        else:
            self.pfc = False
        self._algorithm = alg

    # Checks if input is a single sequence or a file with multiple sequences and creates either a SeqRecord or an iterator
    @property
    def input(self):
        return self._input

    # Input setter function, accepts a Path, string, SeqRecord or List and always returns an Iterator over SeqRecord objects (even if the length is 1)
    @input.setter
    def input(self, input: Path | str | SeqRecord | list[str | SeqRecord]):
        match input:
            case Path():
                if input.is_file():
                    self._input = self._read_input_file(input)
                else:
                    raise LookupError("Could not find input file.")
            case str():
                if Path(input).is_file():
                    self._input = self._read_input_file(input)
                else:
                    if any(c not in "AUCGTaucgt" for c in set(input)):
                        raise ValueError(
                            "Input string was neither a viable file path nor a viable RNA or DNA sequence"
                        )
                    else:
                        self._input = SeqRecord(seq=Seq(input), id="N/A")

            case SeqRecord() | SeqIO.FastaIO.FastaIterator() | SeqIO.QualityIO.FastqPhredIterator():
                self._input = input
            case list():
                if any(type(candidate) is not SeqRecord for candidate in input):
                    raise ValueError("At least one of your list items is not a SeqRecords")

                else:
                    self._input = iter(input)
            case _:
                raise TypeError(
                    "Did not recognize input type, only path, str, SeqRecord, List of Strings or List of SeqRecords is allowed."
                )

    @property
    def call(self):
        """Automatic call setter, if a custom_call is set this will always return the custom_call aswell. The function checks all set runtime arguments in builds a call string"""
        if hasattr(self, "custom_call"):
            return self.custom_call
        runtime_dictionary = {
            "-Q": self.motif_source,
            "-b": self.motif_orientation,
            "-t": self.temperature,
            "-u": self.allowLonelyBasepairs,
        }
        if self.subopt:
            runtime_dictionary["-e"] = self.energy
            runtime_dictionary["-c"] = self.energy_percent
        else:
            runtime_dictionary["-k"] = self.kvalue
        if self.algorithm == "motshapeX":
            runtime_dictionary["-q"] = self.shape_level
        arguments = [
            "{key} {value}".format(key=x, value=y)
            for x, y in runtime_dictionary.items()
            if y is not None
        ]

        call = " ".join([self.algorithm_path, " ".join(arguments), ""])
        if self.time:
            call = "time " + call
        return call

    # Finds File type based on file ending
    def _find_filetype(self, file_path: Path) -> None:
        if str(file_path).split(".")[-1] == "gz" or str(file_path).split(".")[-1] == "zip":
            file_extension = str(file_path).split(".")[-2]
            input_zipped = True
        else:
            file_extension = str(file_path).split(".")[-1]
            input_zipped = False

        match file_extension:
            case "fasta" | "fas" | "fa" | "fna" | "ffn" | "faa" | "mpfa" | "frn" | "txt" | "fsa":
                filetype = "fasta"

            case "fastq" | "fq":
                filetype = "fastq"

            case "stk" | "stockholm" | "sto":
                filetype = "stockholm"
            case _:
                raise TypeError(
                    "Filetype was not recognized as fasta, fastq or stockholm format. Or file could not be unpacked, please ensure it is zipped with either .gz or .zip or unzipped"
                )
        return (input_zipped, filetype)

    # Read input file
    def _read_input_file(
        self, file_path: Path
    ) -> (
        SeqIO.FastaIO.FastaIterator
        | SeqIO.QualityIO.FastqPhredIterator
        | Generator[SeqIO.SeqRecord, None, None]
    ):
        (zipped, filetype) = self._find_filetype(file_path)
        if not zipped:
            return SeqIO.parse(file_path, filetype)
        else:
            with gzip.open(file_path, "rt") as handle:
                return SeqIO.parse(handle, filetype)

    # Calibrate results.algorithm objects based on the current status of the bgap_rna class instance
    def _calibrate_result_objects(self, sep):
        """Calibrate result objects to current configuration (separator, algorithm and pfc + probability filtering)"""
        results.result.separator = sep
        results.result.algorithm = self.algorithm
        results.algorithm_output.pfc = self.pfc
        results.algorithm_output.filtering = self.pfc_filtering

    # Calibrate result objects and run either a single process if the input is a SeqRecord Object or run multiple predictions if input was a file or list of SeqRecord objects
    def auto_run(
        self,
        o_file: Optional[Path] = None,
        pool_workers: int = multiprocessing.cpu_count(),
        output_csv_separator="\t",
    ) -> list[results.algorithm_output]:
        """Checks type of self.input and runs a Single Process in case of a SeqRecord or a MultiProcess in case of an Iterable as input."""
        if output_csv_separator == r"\t":
            output_csv_separator = output_csv_separator.replace(r"\t", "\t")
        self._calibrate_result_objects(output_csv_separator)
        if isinstance(self.input, SeqRecord):
            return self._run_single_process(o_file)
        elif isinstance(
            self.input,
            Iterator
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Generator[SeqIO.SeqRecord, None, None],
        ):
            return self._run_multi_process(o_file, workers=pool_workers)
        else:
            raise TypeError(
                "Input does not match any of the allowed types SeqRecord, FastaIterator, FastqPhredIterator, Iterator[list[SeqRecord]]"
            )

    # single process function utilizing subprocess to run a single prediction and return the output
    def _run_single_process(self, output_f: Path = None) -> results.algorithm_output:
        """Single sequence input running function, silently transcribes DNA to RNA"""
        if "T" in str(self.input.seq):
            self.input.seq = self.input.seq.transcribe()
        subproc_out = subprocess.run(
            self.call + str(self.input.seq),
            text=True,
            capture_output=True,
            shell=True,
            timeout=None,
        )
        if not subproc_out.returncode:
            return_val = results.algorithm_output(
                self.input.id, subproc_out.stdout, subproc_out.stderr
            )
            if isinstance(output_f, Path):
                with open(output_f, "a+") as file:
                    with redirect_stdout(file):
                        return_val.write_results(initiated=False)
            else:
                return_val.write_results(initiated=False)
                sys.stdout.flush()
            return return_val
        else:
            if subproc_out.returncode == 137:
                print(
                    "Your prediction could not be computed due to insufficient memory capacity, please download more RAM."
                )
            else:
                raise subproc_out.stderr

    # multiprocessing function, which utilizes a worker pool to run the specified algorithm on each sequences in the input iterable
    def _run_multi_process(
        self,
        output_file: Optional[Path] = None,
        workers: int = multiprocessing.cpu_count(),
    ) -> Optional[list[results.algorithm_output]]:

        Manager = multiprocessing.Manager()  # Spawn multiprocessing Manager object
        # Start with setting up the input and output Queues
        input_q = Manager.Queue(maxsize=workers * 2)  # type:multiprocessing.Queue[SeqIO.SeqRecord]
        for record in self.input:
            input_q.put(record)

        # Check for the size of the input Q, if it is below workers we have less sequences than workers
        if input_q.qsize() < workers:
            workers = input_q.qsize()

        output_q = Manager.Queue()  # type:multiprocessing.Queue[tuple]
        # Start worker Pool and the listener subprocess, which takes in the worker outputs and formats them
        Pool = multiprocessing.Pool(processes=workers)
        listening = multiprocessing.Process(target=self._listener, args=(output_q, output_file))
        listening.start()
        worker_list = []  # type:list[multiprocessing.AsyncResult]

        for i in range(workers):  # Populate the pool
            work = Pool.apply_async(bgap_rna._worker_funct, (self.call, input_q, output_q))
            worker_list.append(work)

        for i in range(workers):  # Tell the workers to stop
            input_q.put(None)
        for work in worker_list:  # Let workers finish their iterables
            work.get()
        Pool.close()
        Pool.join()
        output_q.put(None)
        listening.join()
        return output_q.get()

    def _listener(
        self, q: multiprocessing.Queue, output_f: Path | None
    ) -> list[results.algorithm_output] | None:
        writing_started = False
        return_list = []
        while True:
            output = q.get()  # type:'results.algorithm_output | results.error | None'
            if output is None:
                q.put(return_list)
                break
            else:
                if isinstance(output, results.algorithm_output):
                    if output_f is not None:
                        with open(output_f, "a+") as file:
                            with redirect_stdout(file):
                                writing_started = output.write_results(writing_started)
                    else:
                        writing_started = output.write_results(writing_started)
                        sys.stdout.flush()
                    return_list.append(output)
                if isinstance(output, results.error):
                    sys.stderr.write(f"{output.id}: {output.error}")
