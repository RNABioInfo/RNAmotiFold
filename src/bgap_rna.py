import sys
import argparse
import subprocess
import configparser
import multiprocessing
from pathlib import Path
from typing import Iterator, Optional, Generator
from Bio import SeqIO
import logging
from Bio.SeqRecord import SeqRecord
from requests import get
import src.results as results
from contextlib import redirect_stdout

import logging

logger = logging.getLogger("bgap_rna")

# A Python class for making Bellman's GAP more convenient to use
# Just create a class instances, feed it with the call arguments you need
# and it'll create a call from the arguments classified as runtime arguments.
# These are: motif_source, motif_orientation, kvalue, hishape_mode, shape_level, energy


class bgap_rna:

    def __repr__(self) -> str:
        return f"bgap_rna({' , '.join([str(x)+'='+str(self.__dict__[x]) for x in self.__dict__])})"

    def __str__(self) -> str:
        return self.call

    @staticmethod
    def _worker_funct(call: str, iq: multiprocessing.Queue, listenerq: multiprocessing.Queue):
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
            listenerq.put(result)

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
        except ConnectionError:
            raise ConnectionError(
                "Could not connect to RNA 3D Motif Atlas server to check for current motif versions."
            )
        if current == installed:
            return False
        else:
            return True

    @classmethod
    def from_argparse(cls, cmd_args: argparse.Namespace):
        """argparse alternative constructor."""
        return cls(
            algorithm=cmd_args.algorithm,
            motif_source=cmd_args.motif_source,
            motif_orientation=cmd_args.motif_orientation,
            kvalue=cmd_args.kvalue,
            shape_level=cmd_args.shape_level,
            energy=cmd_args.energy,
            hishape_mode=cmd_args.hishape_mode,
            temperature=cmd_args.temperature,
            energy_percent=cmd_args.energy_percent,
            allowLonelyBasepairs=cmd_args.basepairs,
            subopt=cmd_args.subopt,
            pfc=cmd_args.pfc,
            pfc_filtering=cmd_args.pfc_filtering,
            session_id=cmd_args.id,
            custom_hairpins=cmd_args.custom_hairpins,
            custom_internals=cmd_args.custom_internals,
            custom_bulges=cmd_args.custom_bulges,
            replace_hairpins=cmd_args.replace_hairpins,
            replace_internals=cmd_args.replace_internals,
            replace_bulges=cmd_args.replace_bulges,
        )

    @classmethod
    def from_config(cls, config: configparser.ConfigParser, section_name: str):
        """ConfigParse alternative constructor. Accepts a string/path to a config file or a ConfigParser. Unused parameters should be left empty or be set to None."""
        return cls(
            algorithm=config.get(section_name, "algorithm"),
            motif_source=config.getint(section_name, "motif_source"),
            motif_orientation=config.getint(section_name, "motif_orientation"),
            kvalue=config.getint(section_name, "kvalue"),
            shape_level=config.getint(section_name, "shape_level"),
            energy=config.get(section_name, "energy"),
            hishape_mode=config.get(section_name, "hishape_mode"),
            temperature=config.getfloat(section_name, "temperature"),
            energy_percent=config.getfloat(section_name, "energy_percent"),
            allowLonelyBasepairs=config.getint(section_name, "basepairs"),
            subopt=config.getboolean(section_name, "subopt"),
            pfc=config.getboolean(section_name, "pfc"),
            pfc_filtering=config.getboolean(section_name, "pfc_filtering"),
            session_id=config.get(section_name, "ID"),
            custom_hairpins=config.get(section_name, "custom_hairpins"),
            custom_internals=config.get(section_name, "custom_internals"),
            custom_bulges=config.get(section_name, "custom_bulges"),
            replace_hairpins=config.getint(section_name, "replace_hairpins"),
            replace_internals=config.getint(section_name, "replace_internals"),
            replace_bulges=config.getint(section_name, "replace_bulges"),
        )

    def __init__(
        self,
        algorithm: Optional[str] = None,
        motif_source: Optional[int] = None,
        motif_orientation: Optional[int] = None,
        kvalue: Optional[int] = None,
        shape_level: Optional[int] = None,
        energy: Optional[str] = None,
        hishape_mode: Optional[str] = None,
        temperature: Optional[float] = None,
        energy_percent: Optional[float] = None,
        custom_hairpins: Optional[str] = None,
        custom_internals: Optional[str] = None,
        custom_bulges: Optional[str] = None,
        allowLonelyBasepairs: bool = False,
        replace_hairpins: bool = False,
        replace_internals: bool = False,
        replace_bulges: bool = False,
        subopt: bool = False,
        pfc: bool = False,
        pfc_filtering: bool = False,
        session_id: str = "N/A",
    ):
        self.id = session_id
        self.subopt = subopt
        self.pfc = pfc
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
        )  # type: SeqIO.FastaIO.FastaIterator | SeqIO.QualityIO.FastqPhredIterator | Generator[SeqRecord, None, None] | Iterator[list] | SeqRecord
        self.pfc_filtering = pfc_filtering  # Set this to enable/disable filtering of partition function outputs to only return results with a probability over 0.0001
        # Custom motif variables, custom_X is for the filepaths to the .csv files, replace_X is for if the customs should append to or replace the underlying motifs from RNA3D or Rfam
        self.custom_hairpins = custom_hairpins
        self.custom_internals = custom_internals
        self.custom_bulges = custom_bulges
        self.replace_hairpins = replace_hairpins
        self.replace_internals = replace_internals
        self.replace_bulges = replace_bulges

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
        if temp is not None:
            if -273 < temp < 100:
                self._temperature = temp
            else:
                raise ValueError("Temperature can only be between -273°C and 100°C")
        else:
            self._temperature = temp

    @property
    def energy_percent(self):
        return self._energy_perc

    @energy_percent.setter
    def energy_percent(self, range: Optional[float]):
        if range is not None:
            if range > 0:
                self._energy_perc = range
            else:
                raise ValueError("Energy range cannot be below 0.")
        elif range is None:
            self._energy_perc = range

    @property
    def allowLonelyBasepairs(self):
        return self._allowLonelyBasepairs

    @allowLonelyBasepairs.setter
    def allowLonelyBasepairs(self, val: bool):
        if val:
            self._allowLonelyBasepairs = 1
        else:
            self._allowLonelyBasepairs = 0

    @property
    def replace_hairpins(self):
        return self._replace_hairpins

    @replace_hairpins.setter
    def replace_hairpins(self, val: bool):
        if val:
            self._replace_hairpins = 1
        else:
            self._replace_hairpins = 0

    @property
    def replace_internals(self):
        return self._replace_internals

    @replace_internals.setter
    def replace_internals(self, val: bool):
        if val:
            self._replace_internals = 1
        else:
            self._replace_internals = 0

    @property
    def replace_bulges(self):
        return self._replace_bulges

    @replace_bulges.setter
    def replace_bulges(self, val: bool):
        if val:
            self._replace_bulges = 1
        else:
            self._replace_bulges = 0

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
                "Motif source can only be 1 = RNA 3D Motif Atlas , 2 = RMfam , 3 = Both"
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
    def energy(self, e: Optional[str]):
        if e == "":
            e = None
        if e is not None:
            if float(e) > 0:
                self._energy = e
            else:
                raise ValueError("Energy range cannot be lower than 0.")
        else:
            self._energy = e

    # Set kvalue fpr kbest and kbacktracing
    @property
    def kvalue(self):
        return self._kvalue

    @kvalue.setter
    def kvalue(self, k: Optional[int]):
        if k is not None:
            if k > 0:
                self._kvalue = k
        elif k is None:
            self._kvalue = k
        else:
            raise ValueError("Kvalue cannot be 0 or lower.")

    # Finds path to your chosen algorithm, if it does not exist i attempts to compile the algorithm
    @property
    def algorithm_path(self):
        return str(
            Path(__file__).resolve().parents[1].joinpath("Build", "bin").joinpath(self.algorithm)
        )

    # Builds the algorithm name string, checks for _pfc at the end and builds hishape name and adds _subopt
    @property
    def algorithm(self):
        return self._algorithm

    @algorithm.setter
    def algorithm(self, alg: Optional[str]):
        if alg is None:
            self._algorithm = None
        else:
            if alg == "RNAmotiCes":
                try:
                    alg = alg + "_" + self.hishape_mode
                except TypeError:
                    raise TypeError("Specify a hishape mode with -p")
            if self.subopt:
                alg = alg + "_subopt"
            elif self.pfc:
                alg = alg + "_pfc"
                if self.subopt:
                    raise ValueError("Partition function can't be used in combination with subopt")
            self._algorithm = alg

    @property
    def call(self):
        """Automatic call setter, if a custom_call is set this will always return the custom_call. The function checks the set algorithm and builds a call string based on it."""
        if hasattr(self, "custom_call"):
            return self.custom_call
        runtime_dictionary = {
            "-Q": self.motif_source,
            "-b": self.motif_orientation,
            "-t": self.temperature,
            "-u": self.allowLonelyBasepairs,
            "-X": self.custom_hairpins,
            "-Y": self.custom_internals,
            "-Z": self.custom_bulges,
            "-L": self.replace_hairpins,
            "-E": self.replace_internals,
            "-G": self.replace_bulges,
        }
        if self.subopt:
            runtime_dictionary["-e"] = self.energy
            runtime_dictionary["-c"] = self.energy_percent
        else:
            runtime_dictionary["-k"] = self.kvalue
        if self.algorithm == "RNAmoSh":
            runtime_dictionary["-q"] = self.shape_level
        arguments = [
            "{key} {value}".format(key=x, value=y)
            for x, y in runtime_dictionary.items()
            if y is not None
        ]

        call = " ".join([self.algorithm_path, " ".join(arguments), ""])
        return call

    @call.setter
    def call(self):
        raise ValueError("Please use the custom_call property so set a custom call.")

    # Calibrate results.algorithm objects based on the current status of the bgap_rna class instance
    def _calibrate_result_objects(
        self, sep: Optional[str] = "\t", logger: Optional[logging.Logger] = None
    ):
        """Calibrate result objects to current configuration (separator, algorithm and pfc + probability filtering)"""
        results.result._separator = sep
        results.result._algorithm = self.algorithm
        results.algorithm_output.pfc = self.pfc
        results.algorithm_output.filtering = self.pfc_filtering

    # Calibrate result objects and run either a single process if the input is a SeqRecord Object or run multiple predictions if input was a file or list of SeqRecord objects
    def auto_run(
        self,
        user_input: (
            SeqRecord
            | Iterator
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Generator[SeqIO.SeqRecord, None, None]
        ),
        o_file: Optional[Path] = None,
        pool_workers: int = multiprocessing.cpu_count(),
        output_csv_separator="\t",
    ) -> list[results.algorithm_output]:
        """Checks type of self.input and runs a Single Process in case of a SeqRecord or a MultiProcess in case of an Iterable as input."""
        if output_csv_separator == r"\t":
            output_csv_separator = output_csv_separator.replace(r"\t", "\t")
        self._calibrate_result_objects(output_csv_separator)
        if isinstance(user_input, SeqRecord):
            return self._run_single_process(user_input, o_file)
        elif isinstance(
            user_input,
            Iterator
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Generator[SeqIO.SeqRecord, None, None],
        ):
            return self._run_multi_process(user_input, o_file, workers=pool_workers)
        else:
            raise TypeError(
                "Input does not match any of the allowed types SeqRecord, FastaIterator, FastqPhredIterator, Iterator[list[SeqRecord]]"
            )

    # single process function utilizing subprocess to run a single prediction and return the output
    def _run_single_process(
        self,
        user_input,
        output_f: Path | None = None,
    ) -> results.algorithm_output | results.error:
        """Single sequence input running function, silently transcribes DNA to RNA"""
        if "T" in str(user_input):
            logger.debug("Detected DNA sequence as input, transcribing to RNA")
            user_input.seq = user_input.seq.transcribe()
            logger.debug(f"Transcription complete, new input:{str(user_input.seq)}")
        logger.debug("Running prediction")
        subproc_out = subprocess.run(
            self.call + str(user_input.seq),
            text=True,
            capture_output=True,
            shell=True,
            timeout=None,
        )
        if not subproc_out.returncode:
            return_val = results.algorithm_output(
                user_input.id, subproc_out.stdout, subproc_out.stderr
            )
            logger.info(f"Prediction finished successfully")
            if isinstance(output_f, Path):
                with open(output_f, "a+") as file:
                    with redirect_stdout(file):
                        return_val.write_results(initiated=False)
                        return return_val
            else:
                return_val.write_results(initiated=False)
                sys.stdout.flush()
            return return_val
        else:
            if subproc_out.returncode == 137:
                print(
                    "Your prediction could not be computed due to insufficient memory capacity on your computer, please download more RAM."
                )
            else:
                logger.info(f"Process {user_input.id} finished with error: {subproc_out.stderr}")
                return results.error(user_input.id, subproc_out.stderr)

    # multiprocessing function, which utilizes a worker pool to run the specified algorithm on each sequences in the input iterable
    def _run_multi_process(
        self,
        user_input: (
            Iterator
            | SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Generator[SeqIO.SeqRecord, None, None]
        ),
        output_file: Optional[Path] = None,
        workers: int = multiprocessing.cpu_count(),
    ) -> Optional[list[results.algorithm_output]]:

        Manager = multiprocessing.Manager()  # Spawn multiprocessing Manager object
        # Start with setting up the input and output Queues
        input_q = Manager.Queue()  # type:multiprocessing.Queue[SeqIO.SeqRecord]
        for record in user_input:
            input_q.put(record)

        # Check for the size of the input Q, if it is below workers we have less sequences than workers
        if input_q.qsize() < workers:
            workers = input_q.qsize()
        listener_q = Manager.Queue()  # type:multiprocessing.Queue[tuple]
        # Start worker Pool and the listener subprocess, which takes in the worker outputs and formats them
        Pool = multiprocessing.Pool(processes=workers)
        listening = multiprocessing.Process(target=self._listener, args=(listener_q, output_file))
        listening.start()
        worker_list = []  # type:list[multiprocessing.AsyncResult]

        for i in range(workers):  # Populate the pool
            work = Pool.apply_async(bgap_rna._worker_funct, (self.call, input_q, listener_q))
            worker_list.append(work)

        for i in range(workers):  # Tell the workers to stop
            input_q.put(None)
        for work in worker_list:  # Let workers finish their iterables
            work.get()
        Pool.close()
        Pool.join()
        listener_q.put(None)
        listening.join()
        return listener_q.get()

    def _listener(
        self,
        q: multiprocessing.Queue,
        output_f: Path | None,
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
                    if isinstance(output_f, Path):
                        logger.info(f"Prediction {output.id} finished.")
                        with open(output_f, "a+") as file:
                            with redirect_stdout(file):
                                writing_started = output.write_results(writing_started)
                    else:
                        writing_started = output.write_results(writing_started)
                        sys.stdout.flush()
                        logger.info(f"Prediction {output.id} finished.")
                    return_list.append(output)
                if isinstance(output, results.error):
                    logger.info(
                        f"Prediction {output.id} finished with error: {output.error.strip()}"
                    )
