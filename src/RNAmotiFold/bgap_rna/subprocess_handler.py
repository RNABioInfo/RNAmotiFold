import multiprocessing
import multiprocessing.connection
import multiprocessing.pool
import subprocess
from typing import Any, Literal
from pathlib import Path
from RNAmotiFold.results.outputs import MFEOutput, PFCOutput, AlignmentOutput, ErrorOutput,AlgorithmOutput
from RNAmotiFold.bgap_rna.input_handler import AlgorithmInput
from RNAmotiFold.input.parameters import ScriptParameters
from contextlib import redirect_stdout
import sys
from collections import defaultdict
import os
import logging

logger = logging.getLogger(__name__)


class SubprocessHandler:
    """Class to manage subprocesses for running RNAmotiFold algorithms."""

    @staticmethod
    def _worker(
        input_queue: "multiprocessing.Queue[AlgorithmInput | None]",
        output_queue: "multiprocessing.Queue[MFEOutput | PFCOutput | AlignmentOutput | ErrorOutput | None]",
        process_type: Literal["mfe", "pfc", "ali"],
    ):
        """Simplest worker function that should work universally, do all pre/post processing outside of this."""
        pid = os.getpid()
        logger.debug(f"Started worker process at {pid}, with process type {process_type}")
        match process_type:
            case "mfe":
                type_alias = MFEOutput
            case "pfc":
                type_alias = PFCOutput
            case "ali":
                type_alias = AlignmentOutput

        while True:
            input_obj: AlgorithmInput | None = input_queue.get()
            if input_obj is None:
                logger.debug(f"Worker {pid} exiting, input queue is empty")
                output_queue.put(None)
                break
            logger.debug(
                f"Worker {pid} started work on {input_obj.id} with call {input_obj.runtime_call}"
            )
            subprocess_output = subprocess.run(
                input_obj.runtime_call, text=True, capture_output=True, shell=True
            )
            logger.debug(f"Worker {pid} processed {input_obj.id}, returned: {subprocess_output}")
            if subprocess_output.returncode == 0:
                result = type_alias(
                    input_obj.id, subprocess_output.stdout, [subprocess_output.stderr]
                )
                logger.debug(
                    f"Worker {pid} successfully completed call {input_obj.call} on {input_obj.id}"
                )
            else:
                result = ErrorOutput(
                    input_obj.id, subprocess_output.stderr, [str(subprocess_output.returncode)]
                )
                logger.debug(
                    f"Worker {pid} encountered exit status {subprocess_output.returncode} working on {input_obj.id}. Error: {subprocess_output.stderr}"
                )
            output_queue.put(result)

    @staticmethod
    def _listener(
        input_queue: "multiprocessing.Queue[MFEOutput | PFCOutput | AlignmentOutput | ErrorOutput | None]",
        output_file: Path | None,
        pipe: multiprocessing.connection.Connection,
        calls_per_input: int,
        no_print: bool,
        max_nones: int,
    ):
        """Simplest listener funtion that should also work universally, takes the result objects put into its queue by the workers, writes them down or prints them.
        When all workers are done, signaled by the sentinel None in the Queue which comes from the main process, terminates and sends a list of result objects to back.
        """
        none_counter = 0
        output_dict: defaultdict[str, list[MFEOutput | PFCOutput | AlignmentOutput]] = defaultdict(
            list
        )  # It will always only be one of the types
        error_dict: defaultdict[str, list[ErrorOutput]] = defaultdict(list)
        writing_started = False
        while True:
            try:
                result: None | MFEOutput | PFCOutput | AlignmentOutput | ErrorOutput = (
                    input_queue.get()
                )
            except EOFError:
                continue
            if result is None:
                none_counter += 1
                if none_counter == max_nones:
                    # Everything is done, pack up all the results and send them back
                    resultlist: list[MFEOutput | PFCOutput | AlignmentOutput] = [
                        item for values in output_dict.values() for item in values
                    ]
                    errorlist: list[ErrorOutput] = [
                        item for values in error_dict.values() for item in values
                    ]
                    pipe.send(resultlist + errorlist)
                    logger.debug("All algorithm returns have been processed, listener is exiting")
                    break
                continue
            match result:
                case MFEOutput() | PFCOutput() | AlignmentOutput():
                    output_dict[result.id].append(result)
                    if (
                        sum([len(output_dict[result.id]), len(error_dict[result.id])])
                        == calls_per_input
                    ):
                        successes = [x for x in output_dict[result.id]]
                        if all([isinstance(x, MFEOutput) for x in successes]):
                            merged = [MFEOutput.merge(successes)]  # type: ignore All 3 of these are checked beforehand and therefore will never fail
                        elif all([isinstance(x, PFCOutput) for x in successes]):
                            merged = PFCOutput.merge(successes)  # type: ignore All 3 of these are checked beforehand and therefore will never fail
                        elif all([isinstance(x, AlignmentOutput) for x in output_dict[result.id]]):
                            merged = [AlignmentOutput.merge(successes)]  # type: ignore All 3 of these are checked beforehand and therefore will never fail
                        else:
                            raise ValueError(
                                "Multiple algorithm output types detected during the same run, something went wrong"
                            )
                        if not no_print:
                            if isinstance(output_file, Path):
                                write_destination = open(output_file, "a+")
                            else:
                                write_destination = sys.stdout
                            with redirect_stdout(write_destination):
                                for element in merged:
                                    writing_started = element.write_results(writing_started)
                                sys.stdout.flush()
                case ErrorOutput():
                    result.write_results(False)
                    error_dict[result.id].append(result)

    def __init__(
        self,
        max_processes: int,
        output_path: Path | None,
        prediction_type: Literal["mfe", "pfc", "ali"],
    ):
        # Set up the very basics of a new multiprocessing step, How Many processes can we start, what are our inputs and where are we supposed to write the output.
        # If we come from the full RNAmotiFold pipeline the input is pre-processed and output location is confirmed to work -> Alternate contsructor for checking these ?
        self.max_processes: int = max_processes
        self.output_path: Path | None = output_path
        self.process_type: Literal["mfe"] | Literal["pfc"] | Literal["ali"] = prediction_type

    @classmethod
    def from_script_parameters(cls, params: ScriptParameters):
        return cls(params.workers, params.output, params.alg_type())

    def auto_run(
        self,
        inputs: list[AlgorithmInput],
        calls_per_input: int,
        no_print: bool,
    ) -> AlgorithmOutput:
        if len(inputs) >= 2 and self.max_processes > 1:
            logger.debug("Recognized multiple inputs, creating worker pool and running in parallel")
            return AlgorithmOutput.from_list(self._multi_run(inputs, calls_per_input, no_print))
        else:
            logger.debug("Recognized only one input, running in single process mode")
            answer = self._single_run(inputs[0])
            if not no_print:
                answer.write_results(False)
            return AlgorithmOutput.from_list([answer])

    def _single_run(
        self,
        input_obj: AlgorithmInput,
    ) -> MFEOutput | PFCOutput | AlignmentOutput | ErrorOutput:
        """Function for running a single"""
        subprocess_output = subprocess.run(
            input_obj.runtime_call, text=True, capture_output=True, shell=True
        )
        if subprocess_output.returncode == 0:
            match self.process_type:
                case "mfe":
                    return MFEOutput(
                        input_obj.id, subprocess_output.stdout, [subprocess_output.stderr]
                    )
                case "pfc":
                    return PFCOutput(
                        input_obj.id, subprocess_output.stdout, [subprocess_output.stderr]
                    )
                case "ali":
                    return AlignmentOutput(
                        input_obj.id, subprocess_output.stdout, [subprocess_output.stderr]
                    )
        else:
            return ErrorOutput(
                input_obj.id, subprocess_output.stderr, [str(subprocess_output.returncode)]
            )

    def _multi_run(
        self,
        inputs: list[AlgorithmInput],
        calls_per_input: int,
        no_print: bool,
    ) -> list[MFEOutput | PFCOutput | AlignmentOutput | ErrorOutput]:
        # Set Up Everything for a multiprocessed run, first make a multiprocessing manager and fill the worker queue with inputs
        manager = multiprocessing.Manager()
        input_q: "multiprocessing.Queue[AlgorithmInput | None]" = manager.Queue()  # type: ignore Because Queue has type Any
        listener_q = manager.Queue()
        if len(inputs) < self.max_processes:
            logger.debug(
                f"Number of inputs is less than max allowed processes ({self.max_processes}), starting only {len(inputs)} workers."
            )
            necessary_processes = len(inputs)
        else:
            necessary_processes = self.max_processes
        for record in inputs:
            input_q.put(record)
        # First we set up  a listener with a connection to the main process
        PipeOut, PipeIn = multiprocessing.Pipe(duplex=False)
        listening = multiprocessing.Process(
            target=self._listener,
            args=(
                listener_q,
                self.output_path,
                PipeIn,
                calls_per_input,
                no_print,
                necessary_processes,
            ),
        )
        listening.start()

        # Now we make a pool of workers and
        pool = multiprocessing.Pool(necessary_processes)
        workers: list[multiprocessing.pool.AsyncResult[Any]] = []
        for _ in range(
            necessary_processes
        ):  # Populate the pool with worker functions, each doing nothing but getting items from the input queue and processing them
            work = pool.apply_async(
                SubprocessHandler._worker, (input_q, listener_q, self.process_type)
            )
            workers.append(work)
            input_q.put(None)

        # Close the Pool
        pool.close()
        pool.join()
        # Tell the listener that no more objects will be put into the queue with the None sentinel and let him finish working
        listening.join()
        listener_output: list[MFEOutput | PFCOutput | AlignmentOutput | ErrorOutput] = (
            PipeOut.recv()
        )  # Receive the list of outputs from the listener
        return listener_output
