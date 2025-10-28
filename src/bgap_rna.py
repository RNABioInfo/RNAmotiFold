import sys
import subprocess
import multiprocessing
from pathlib import Path
from typing import Optional, Generator, Any, Literal 
from Bio.SeqIO import FastaIO, QualityIO
import logging
from Bio.SeqRecord import SeqRecord
from src.args import script_parameters
from src.results import algorithm_output,error
import src.results as results
from contextlib import redirect_stdout

logger = logging.getLogger("bgap_rna")

# A Python class for making Bellman's GAP more convenient to use
# Just create a class instances, feed it with the call arguments you need
# and it'll create a call from the arguments classified as runtime arguments.
# These are: motif_source, motif_orientation, kvalue, hishape_mode, shape_level, energy


class bgap_rna:
    '''
    Main class for running RNAmotiFold algortihms through python, just hand your arguments to this class (everything else will be defaults) and use [your class obj].auto_run([input]) to run predictions.)
    To be agile with 
    '''
    def __repr__(self):
        classname = type(self).__name__
        k, v = zip(*self.__dict__.items())
        together: list[str] = []
        for i in range(0, len(v)):
            together.append("{key}={value!r}".format(key=k[i], value=v[i]))
        return f"{classname}({', '.join(together)})"

    def __str__(self) -> str:
        return self.call

    @staticmethod
    def _worker(
        call: str,
        iq: "multiprocessing.Queue[SeqRecord|None]",
        listenerq: "multiprocessing.Queue[results.algorithm_output|results.error]",
    ):
        while True:
            record: Optional[SeqRecord] = iq.get()
            if record is None:
                break
            else:
                if "T" in str(record.seq):
                    record.seq = record.seq.transcribe() #type:ignore This throws a transcribe is not Part of none error, though record.seq can't be None at this point 
                logger.debug(f"Running prediction: {record.id}:{call+str(record.seq)}")
                subprocess_output = subprocess.run(
                    call + str(record.seq), text=True, capture_output=True, shell=True
                )
            if not subprocess_output.returncode:
                result = results.algorithm_output(
                    str(record.id), subprocess_output.stdout, [subprocess_output.stderr]
                )
            else:
                result = results.error(str(record.id), subprocess_output.stderr)
            listenerq.put(result)

    @staticmethod
    def _listener(
        q: "multiprocessing.Queue[results.algorithm_output | results.error | list[results.algorithm_output|results.error] | None]",
        output_f: Optional[Path | str] = None,
    ) -> None:
        writing_started = False
        return_list: list[results.algorithm_output | results.error] = []
        while True:
            output: None | results.algorithm_output | results.error | list[results.algorithm_output | results.error]= q.get()
            if output is None:
                q.put(return_list)
                break
            else:
                if isinstance(output, results.algorithm_output):
                    logger.info(f"Prediction {output.id} finished successfully.")
                    if isinstance(output_f, (Path,str)):
                        with open(output_f, "a+") as file:
                            with redirect_stdout(file):
                                writing_started = output.write_results(writing_started)
                    else:
                        writing_started = output.write_results(writing_started)
                        sys.stdout.flush()   
                    return_list.append(output)
                if isinstance(output, results.error):
                    logger.warning(
                        f"Prediction {output.id} finished with error: {output.error.strip()}"
                    )
                    return_list.append(output)

    @classmethod
    def from_script_parameters(cls, params: script_parameters):
        return cls(
            alg=params.algorithm,
            motif_source=params.motif_source,
            motif_orientation=params.motif_orientation,
            kvalue=params.kvalue,
            shape_level=params.shape_level,
            energy=params.energy,
            pfc=params.pfc,
            low_prob_filter=params.low_prob_filter,
            temperature=params.temperature,
            energy_percent=params.energy_percent,
            custom_hairpins=params.custom_hairpins,
            custom_internals=params.custom_internals,
            custom_bulges=params.custom_bulges,
            allowLonelyBasepairs=params.basepairs,
            replace_hairpins=params.replace_hairpins,
            replace_internals=params.replace_internals,
            replace_bulges=params.replace_bulges,
            subopt=params.subopt,
            session_id=params.id,
            fast_mode=params.fast_mode,
        )

    def __init__(
        self,
        alg: Literal["RNAmotiFold","RNAmoSh","RNAmotiCes"],
        motif_source: int = 1,
        motif_orientation: int = 1,
        kvalue: int = 5,
        shape_level: int = 3,
        energy: Optional[str] = None,
        temperature: float = 37.0,
        energy_percent: float = 5.0,
        custom_hairpins: Optional[Path|str] = None,
        custom_internals: Optional[Path|str] = None,
        custom_bulges: Optional[Path|str] = None,
        allowLonelyBasepairs: bool = False,
        replace_hairpins: Optional[bool] = None,
        replace_internals: Optional[bool] = None,
        replace_bulges: Optional[bool] = None,
        subopt: bool = False,
        pfc: bool = False,
        low_prob_filter: float = 0.000001,
        session_id: str = "N/A",
        fast_mode:bool = False,
    ):
        self.id = session_id
        self.subopt = subopt
        self.pfc = pfc
        self.low_probability_filter = low_prob_filter
        self.motif_source = motif_source
        self.motif_orientation = motif_orientation
        self.kvalue = kvalue
        self.shape_level = shape_level
        self.absolute_energy = energy
        self.temperature = temperature
        self.algorithm = alg
        self.energy_percent = energy_percent
        self.allowLonelyBasepairs = allowLonelyBasepairs
        # Custom motif variables, custom_X is for the filepaths to the .csv files, replace_X is for if the customs should append to or replace the underlying motifs from RNA3D or Rfam
        self.custom_hairpins = custom_hairpins
        self.custom_internals = custom_internals
        self.custom_bulges = custom_bulges
        self.replace_hairpins = replace_hairpins
        self.replace_internals = replace_internals
        self.replace_bulges = replace_bulges
        self.fast_mode = fast_mode

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
    def temperature(self, temp: float):
        if -273 < temp < 100:
            self._temperature = temp
        else:
            raise ValueError("Temperature can only be between -273°C and 100°C")

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
    def replace_hairpins(self, val: Optional[bool]):
        if val is None:
            self._replace_hairpins = None
        elif val:
            self._replace_hairpins = 1
        else:
            self._replace_hairpins = 0

    @property
    def replace_internals(self):
        return self._replace_internals

    @replace_internals.setter
    def replace_internals(self, val: Optional[bool]):
        if val is None:
            self._replace_internals = None
        elif val:
            self._replace_internals = 1
        else:
            self._replace_internals = 0

    @property
    def replace_bulges(self):
        return self._replace_bulges

    @replace_bulges.setter
    def replace_bulges(self, val: Optional[bool]):
        if val is None:
            self._replace_bulges = None
        elif val:
            self._replace_bulges = 1
        else:
            self._replace_bulges = 0

    # Choose motif source, 1 = RNA 3D Motif Atlas, 2 = RMFam, 3 = Both
    @property
    def motif_source(self):
        return self._motif_source

    @motif_source.setter
    def motif_source(self, Q: int):
        if Q in [1, 2, 3]:
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
    def motif_orientation(self, b: int):
        if b in [1, 2, 3]:
            self._motif_orientation = b
        else:
            raise ValueError("Motif direction can only be 1 = 5'->3' , 2 = 3'->5' , 3 = Both.")

    # Set shape abstraction level, viable inputs are 1-5
    @property
    def shape_level(self):
        return self._shape_level

    @shape_level.setter
    def shape_level(self, q: int):
        if q in [1, 2, 3, 4, 5]:
            self._shape_level = q
        else:
            raise ValueError(
                "Shape level can only be set to levels 5 (most abstract) - 1 (least abstract)."
            )

    # Set energy range for suboptimal candidate calcualtions
    @property
    def absolute_energy(self):
        return self._energy

    @absolute_energy.setter
    def absolute_energy(self, e: Optional[str]):
        if e == "":
            e = None
        if e is not None:
            if float(e) >= 0:
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
    def kvalue(self, k: int):
        if k > 0:
            self._kvalue = k
        else:
            raise ValueError("Kvalue cannot be 0 or lower.")

    @property
    def low_probability_filter(self):
        return self._low_probability_filter

    @low_probability_filter.setter
    def low_probability_filter(self, value: float):
        if 0 < value < 1:
            self._low_probability_filter = value
        else:
            raise ValueError("Probability filter cannot be below 0 or above 1")

    # Finds path to your chosen algorithm, if it does not exist i attempts to compile the algorithm
    @property
    def algorithm_path(self):
        return str(
            Path(__file__).resolve().parents[1].joinpath("Build", "bin").joinpath(self.algorithm)
        )

    # Builds the algorithm name string, adds _pfc or _subopt
    @property
    def algorithm(self):
        return self._algorithm

    @algorithm.setter
    def algorithm(self, alg: str):
        if self.subopt:
            alg = alg + "_subopt"
        elif self.pfc:
            alg = alg + "_pfc"
            if self.subopt:
                raise ValueError("Partition function can't be used in combination with subopt")
        self._algorithm = alg

    @property
    def custom_hairpins(self) -> Path|None:
        return self._custom_hairpins
    
    @custom_hairpins.setter
    def custom_hairpins(self,path:str|Path|None) -> None:
        if path is not None:
            if Path(path).is_file():
                self._custom_hairpins = Path(path)
            else:
                raise FileNotFoundError("Unabled to find specified custom hairpin file, please check the file path")
        else:
            self._custom_hairpins = None

    @property
    def custom_internals(self) -> Path|None:
        return self._custom_internals
    
    @custom_internals.setter
    def custom_internals(self,path:str|Path|None) -> None:
        if path is not None:
            if Path(path).is_file():
                self._custom_internals = Path(path)
            else:
                raise FileNotFoundError("Unabled to find specified custom internals file, please check the file path")
        else:
            self._custom_internals = None

    @property
    def custom_bulges(self) -> Path|None:
        return self._custom_bulges
    
    @custom_bulges.setter
    def custom_bulges(self,path:str|Path|None) -> None:
        if path is not None:
            if Path(path).is_file():
                self._custom_bulges = Path(path)
            else:
                raise FileNotFoundError("Unabled to find specified custom bulges file, please check the file path")
        else:
            self._custom_bulges = None

    @property
    def call(self):
        """Automatic call setter, if a custom_call is set this will always return the custom_call. The function checks the set algorithm and builds a call string based on it."""
        if hasattr(self, "custom_call"):
            return self.custom_call
        runtime_dictionary: dict[str, Optional[str] | Optional[int] | Optional[float] | Path] = {
            "-F": self.low_probability_filter,
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
        }  # type:dict[str,Optional[float,int,str]]
        if self.subopt:
            # Ordering here is important, the last one is always used so to keep -e overwriting -c this is necessary
            runtime_dictionary["-c"] = self.energy_percent
            runtime_dictionary["-e"] = self.absolute_energy
        else:
            runtime_dictionary["-k"] = self.kvalue
        if self.algorithm in ["RNAmoSh", "RNAmoSh_subopt", "RNAmoSh_pfc"]:
            runtime_dictionary["-q"] = self.shape_level
        arguments = [
            "{key} {value}".format(
                key=x,
                value=y,
            )
            for x, y in runtime_dictionary.items()
            if y is not None
        ]

        seq_free_call = " ".join(
            [self.algorithm_path, " ".join(arguments), ""] #"/usr/bin/time"
        )  # Creates call string without a sequence
        return seq_free_call

    @call.setter
    def call(self):
        raise ValueError("Please use the custom_call property so set a custom call.")

    # Calibrate results.algorithm objects based on the current status of the bgap_rna class instance
    def _calibrate_result_objects(self, sep: str = "\t"):
        """Calibrate result objects to current configuration (separator, algorithm and pfc + probability filtering)"""
        results.result.separator = sep
        if self.pfc:
            results.algorithm_output.Status = "pfc"
        else:
            results.algorithm_output.Status = "mfe"

    # Calibrate result objects and run either a single process if the input is a SeqRecord Object or run multiple predictions if input was a file or list of SeqRecord objects
    def auto_run(
        self,
        user_input: (
            SeqRecord
            | list[SeqRecord]
            | FastaIO.FastaIterator
            | QualityIO.FastqPhredIterator
            | Generator[SeqRecord, None, None]
        ),
        version:str,
        o_file: Optional[Path | str] = None,
        pool_workers: int = multiprocessing.cpu_count(),
        output_csv_separator: str = "\t",
        merge:bool = False,
    ) -> list[results.algorithm_output | results.error ]:
        """Checks type of self.input and runs a Single Process in case of a SeqRecord or a MultiProcess in case of an Iterable as input."""

        if output_csv_separator == r"\t":
            output_csv_separator = output_csv_separator.replace(r"\t", "\t")
        self._calibrate_result_objects(output_csv_separator)
        if self.fast_mode:
            return self.run_separate_processes(user_input,version,o_file,pool_workers,merge)
        if isinstance(user_input, SeqRecord):
            return self._run_single_process(user_input, o_file)
        else:
            return self._run_multi_process(user_input, o_file, workers=pool_workers)

    # single process function utilizing subprocess to run a single prediction and return the output
    def _run_single_process(
        self,
        user_input: SeqRecord,
        output_f: Optional[Path|str] = None,
    ) -> list[results.algorithm_output | results.error]:
        """Single sequence input running function, silently transcribes DNA to RNA"""
        if "T" in str(user_input):
            logger.info("Detected DNA sequence as input, transcribing to RNA")
            user_input.seq = user_input.seq.transcribe() #type:ignore Same case as above, throws transcribe is not part of none error but it can't be None here
            logger.info(f"Transcription complete, new input:{str(user_input.seq)}")
        logger.debug(f"Running prediction:{self.id}:{self.call} {str(user_input.seq)}")
        subproc_out = subprocess.run(
            self.call + str(user_input.seq),
            text=True,
            capture_output=True,
            shell=True,
            timeout=None,
        )
        if not subproc_out.returncode:
            return_val = results.algorithm_output(
                str(user_input.id), subproc_out.stdout, [subproc_out.stderr]
            )
            logger.info(f"Prediction finished successfully")
            if isinstance(output_f, Path):
                with open(output_f, "a+") as file:
                    with redirect_stdout(file):
                        return_val.write_results(initiated=False)
                        return [return_val]
            else:
                return_val.write_results(initiated=False)
                sys.stdout.flush()
            return [return_val]
        else:
            logger.warning(f"Process {user_input.id} finished with error: {subproc_out.stderr}")
            return [results.error(str(user_input.id), subproc_out.stderr)]

    # multiprocessing function, which utilizes a worker pool to run the specified algorithm on each sequences in the input iterable
    def _run_multi_process(
        self,
        user_input: (
            list[SeqRecord]
            | FastaIO.FastaIterator
            | QualityIO.FastqPhredIterator
            | Generator[SeqRecord, None, None]
        ),
        output_file: Optional[Path| str] = None,
        workers: int = multiprocessing.cpu_count(),
    ) -> list[results.algorithm_output | results.error]:

        Manager = multiprocessing.Manager()  # Spawn multiprocessing Manager object
        # Start with setting up the input and output Queues
        input_q: multiprocessing.Queue[SeqRecord | str | None] = Manager.Queue()  # type: ignore
        for record in user_input:
            input_q.put(record)

        # Check for the size of the input Q, if it is below workers we have less sequences than workers
        if input_q.qsize() < workers:
            workers = input_q.qsize()
        listener_q: multiprocessing.Queue[results.algorithm_output | results.error | None | list[results.algorithm_output | results.error]] = Manager.Queue()  # type: ignore
        # Start worker Pool and the listener subprocess, which takes in the worker outputs and formats them
        Pool = multiprocessing.Pool(processes=workers)
        listening = multiprocessing.Process(target=self._listener, args=(listener_q, output_file))
        listening.start()
        worker_list: list[Any] = []

        for i in range(workers):  #type:ignore Populate the pool
            work = Pool.apply_async(bgap_rna._worker, (self.call, input_q, listener_q))
            worker_list.append(work)
            input_q.put(None)
    
        for work in worker_list:  # Let workers finish
            work.get()
        Pool.close()
        Pool.join()
        listener_q.put(None)
        listening.join()
        listener_results: list[results.algorithm_output | results.error] = listener_q.get() #type:ignore This will always be a list of results since it's what the listener "returns"
        return listener_results
    
    @staticmethod
    def get_motif_files(version:str)->list[Path]:
        motif_dir_path = Path(__file__).resolve().parent.joinpath("..","submodules","RNALoops","Misc","Applications","RNAmotiFold","motifs","versions",f"{version}_separated")
        files = Path(motif_dir_path).rglob("*.csv")
        return list(files)
    
    @staticmethod
    def set_customs(bgap_obj:'bgap_rna',motif_path:Path) -> tuple['bgap_rna', Literal['hairpin', 'internal', 'bulge']]:
        path_to_empty_file=Path(__file__).resolve().parent.joinpath("..","submodules","RNALoops","Misc","Applications","RNAmotiFold","motifs","versions","empty.csv")
        match motif_path.parent.name:
            case "hairpins":
                bgap_obj.custom_hairpins = motif_path
                bgap_obj.custom_internals = path_to_empty_file
                bgap_obj.custom_bulges = path_to_empty_file
                motif_type ="hairpin"
            case "internals":
                bgap_obj.custom_hairpins = path_to_empty_file
                bgap_obj.custom_internals = motif_path
                bgap_obj.custom_bulges = path_to_empty_file
                motif_type="internal"
            case "bulges":
                bgap_obj.custom_hairpins = path_to_empty_file
                bgap_obj.custom_internals = path_to_empty_file
                bgap_obj.custom_bulges = motif_path
                motif_type="bulge"
            case _:
                raise FileNotFoundError("Something went wrong with the files paths, unable to find motif files")
        return (bgap_obj, motif_type)
    
    @staticmethod
    def sep_motif_worker_func(bgap_obj:'bgap_rna',input_q:'multiprocessing.JoinableQueue[tuple[SeqRecord,Path]|None]',listener_q:'multiprocessing.Queue[results.algorithm_output|results.error]') -> None:
        #Initialization of each worker, set all replacements to True, assign the right path
        bgap_obj.replace_hairpins = True
        bgap_obj.replace_internals = True
        bgap_obj.replace_bulges = True
        while True:
            record= input_q.get()
            if record is None:
                input_q.task_done()
                break
            record,path = record
            (bgap_obj,motif_type) = bgap_rna.set_customs(bgap_obj,path)
            record.seq = record.seq.transcribe() #type:ignore
            subprocess_output = subprocess.run(bgap_obj.call + str(record.seq),text=True,capture_output=True,shell=True)
            if not subprocess_output.returncode:
                res = results.algorithm_output(str(record.id),subprocess_output.stdout,[subprocess_output.stderr],motif_type)
            else:
                res = results.error(str(record.id),subprocess_output.stderr)
            listener_q.put(res)
            input_q.task_done()

    @staticmethod
    def postprocessing(merged_output:results.algorithm_output) -> results.algorithm_output:
        '''
        Postprocessing function for merging outputs of the seperated motif predictions
        '''
        mfe_dict:dict[int,list[results.result_mfe]] = {}
        for res in merged_output.results:
            if isinstance(res,results.result_mfe) and res.classifier: #this is a little unnecessary but it gets rid of warnings, the res classifier filter removes the "no motif" structure
                if res.free_energy not in mfe_dict.keys():                               #-> It makes no sense to have it in the merging process since if it can fit a motif it will be the mfe for that motif anyways
                    mfe_dict[res.free_energy] = [res]
                else:
                    mfe_dict[res.free_energy].append(res)
        for key in mfe_dict.keys():
            if len(mfe_dict[key]) > 1:
                merge_candidates = results.result_mfe.get_compatible_structures(mfe_dict[key])
                for compatible_structures in merge_candidates:
                    new_result = results.result_mfe.merge_structures([mfe_dict[key][i] for i in compatible_structures])
                    if new_result is not None:
                        merged_output.results.append(new_result)
            else:
                continue
        merged_output.results.sort(key=lambda x: int(x.free_energy)) #type:ignore
        return merged_output

    @staticmethod
    def _listener_sep_motifs(motif_amount:int,q: 'multiprocessing.Queue[results.algorithm_output | results.error | list[results.algorithm_output|results.error] | None]',output_f: Path | None,merge:bool = False) -> None:
        output_dict:dict[str,list[results.algorithm_output]] = {}
        writing_started = False
        return_list: list[results.algorithm_output|results.error] = []
        while True:
            output: 'results.algorithm_output | results.error | list[results.algorithm_output | results.error] | None' = q.get()
            if output is None:
                q.put(return_list)
                break
            else:
                if isinstance(output, results.algorithm_output):
                    if output.id not in output_dict.keys():
                        output_dict[output.id] = [output]
                    else:
                        output_dict[output.id].append(output)
                    if len(output_dict[output.id]) == motif_amount:
                        logger.info(f"Motif predictions for sequence {output.id} finished, merging results.")
                        output = results.algorithm_output.merge_outputs(output_dict[output.id])
                        if merge:
                            logger.info(f"Merging finished for {output.id}, postprocessing...")
                            output = bgap_rna.postprocessing(output)
                        if isinstance(output_f, Path):
                            with open(output_f, "a+") as file:
                                with redirect_stdout(file):
                                    writing_started = output.write_results(writing_started)
                        else:
                            writing_started = output.write_results(writing_started)
                            sys.stdout.flush()
                        return_list.append(output)
                    if isinstance(output, results.error):
                        return_list.append(output)

    def run_separate_processes(self,user_input: (list[SeqRecord]| FastaIO.FastaIterator| QualityIO.FastqPhredIterator| Generator[SeqRecord, None, None] | SeqRecord), version:str, output_file: Optional[Path| str] = None,workers: int = multiprocessing.cpu_count(),merge:bool = False) -> list[algorithm_output | error]:
        '''
        Fast mode multi processing function, fast mode is always multi-processed since parallelizing motifs is what makes it fast (duh).
        '''
        #First we get the motif files we need from our motif sets (I should also change this in get_motif_sequences down the line, concat the files and make a hexdump that way for normal)
        file_list = bgap_rna.get_motif_files(version=version)

        #setup shared object manager for communication between processes, with a list of input sequences
        Manager = multiprocessing.Manager()
        input_q: multiprocessing.JoinableQueue[tuple[SeqRecord,Path]|None] = Manager.JoinableQueue() #type:ignore

        #setup listener and listener queue --> put postprocessing into listener
        listener_q: multiprocessing.Queue[results.algorithm_output | results.error | None | list[results.algorithm_output | results.error]] = Manager.Queue()  # type: ignore Weird Any typing of Queue throws type error
        listening = multiprocessing.Process(target=self._listener_sep_motifs, args=(len(file_list),listener_q, output_file,merge))
        listening.start()

        #Populate the input_list, allocate file paths to each sequences so every sequence gets run with every motif file
        if not isinstance(user_input,SeqRecord):
            for seq in user_input:
                for file in file_list:
                    pair = seq,file
                    input_q.put(pair)
        else:
            for file in file_list:
                pair = user_input,file
                input_q.put(pair)

        #Populate pool with the workers
        Pool = multiprocessing.Pool(processes=workers)
        worker_list2:list[Any] = []

        for _ in range(workers):
            work = Pool.apply_async(func=bgap_rna.sep_motif_worker_func,args=(self,input_q,listener_q))
            worker_list2.append(work)

        for _ in range(len(worker_list2)):
            input_q.put(None)

        #Wait for all tasks to be done and close everything
        input_q.join()
        Pool.close()
        Pool.join()
        listener_q.put(None) #Shuts down the listener after all work is done
        listening.join()
        listener_results: list[results.algorithm_output|results.error] = listener_q.get() #type:ignore This will always be a list of results since the listener puts a list of results into the queue when it shuts down
        return listener_results