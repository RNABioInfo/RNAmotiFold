from __future__ import annotations

import logging
import multiprocessing
import multiprocessing.connection
import multiprocessing.pool
import os
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, NamedTuple, Self

from RNAmotiFold.bgap_rna.input_handler import (
    ReadyAlgorithmInput,
)
from RNAmotiFold.input.parameters import ScriptParameters
from RNAmotiFold.Misc.EnumTypes import PredictionType
from RNAmotiFold.Misc.Helpers import flatten
from RNAmotiFold.results.outputs import (
    AlgorithmOutput,
    AlignmentOutput,
    ErrorOutput,
    MFEOutput,
    PFCOutput,
)

logger = logging.getLogger(__name__)


def create_subprocess_handler(
    params: ScriptParameters,
) -> SubprocessHandler:
    return make_handler(params.workers, params.output, params.alg_type)


def make_handler(
    max_procs: int,
    output: Path | None,
    algorithm_type: PredictionType,
) -> SubprocessHandler:
    match algorithm_type:
        case PredictionType.MFE:
            return SubprocessHandler(max_procs, output, MFEOutput)
        case PredictionType.PFC:
            return SubprocessHandler(max_procs, output, PFCOutput)
        case PredictionType.ALI:
            return SubprocessHandler(max_procs, output, AlignmentOutput)


class CallandID(NamedTuple):
    call: list[str]
    id: str


class SubprocessHandler[T: (MFEOutput, PFCOutput, AlignmentOutput)]:
    """Class to manage subprocesses for running RNAmotiFold algorithms."""

    def __init__(
        self,
        max_processes: int,
        output_path: Path | None,
        result_type: type[T],
    ) -> None:
        # Set up the very basics of a new multiprocessing step, How Many processes can we start, what are our inputs and where are we supposed to write the output.
        # If we come from the full RNAmotiFold pipeline the input is pre-processed and output location is confirmed to work -> Alternate contsructor for checking these ?
        self.max_processes: int = max_processes
        self.output_path: Path | None = output_path
        self.result_type = result_type

    @classmethod
    def from_script_parameters(cls, params: ScriptParameters, type: type[T]) -> Self:
        return cls(params.workers, params.output, type)

    @staticmethod
    def _worker(
        input_queue: multiprocessing.Queue[CallandID | None],
        output_queue: multiprocessing.Queue[T | ErrorOutput | None],
        return_type: type[T],
    ) -> None:
        """Simplest worker function that should work universally, do all pre/post processing outside of this."""
        pid = os.getpid()
        logger.debug(f"Started worker process at {pid}")
        while True:
            ReadyInput: CallandID | None = input_queue.get()
            if ReadyInput is None:
                logger.debug(f"Worker {pid} exiting, input queue is empty")
                output_queue.put(None)
                break
            logger.debug(
                f"Worker {pid} started work on {ReadyInput.id} with call {' '.join(ReadyInput.call)}"
            )
            try:
                subprocess_output = subprocess.run(
                    ReadyInput.call,
                    text=True,
                    capture_output=True,
                    check=True,
                )
                logger.debug(
                    f"Worker {pid} successfully completed call {ReadyInput.call} on {ReadyInput.id}"
                )
                result = return_type(
                    ReadyInput.id,
                    subprocess_output.stdout,
                    [subprocess_output.stderr],
                )
            except subprocess.CalledProcessError as exc:
                result = ErrorOutput(
                    ReadyInput.id,
                    " ".join(ReadyInput.call),
                    exc.stderr,
                    [str(exc.returncode)],
                )
                logger.debug(
                    f"Worker {pid} encountered exit status {exc.returncode} working on {ReadyInput.id}. Error: {exc.stderr}"
                )
            except OSError as oserr:
                result = ErrorOutput(
                    ReadyInput.id,
                    " ".join(ReadyInput.call),
                    str(oserr.strerror),
                    [str(oserr.errno)],
                )
            output_queue.put(result)

    @staticmethod
    def _listener(
        input_queue: multiprocessing.Queue[T | ErrorOutput | None],
        output_file_path: Path | None,
        pipe: multiprocessing.connection.Connection,
        calls_per_input: int,
        no_print: bool,
        max_nones: int,
        typer: type[T],
    ) -> None:
        """Simplest listener funtion that should also work universally, takes the result objects put into its queue by the workers, writes them down or prints them.
        When all workers are done, signaled by the sentinel None in the Queue which comes from the main process, terminates and sends a list of result objects to back.
        """
        none_counter = 0
        output_dict: defaultdict[str, list[T | ErrorOutput]] = defaultdict(list)
        writing_started = False
        while True:
            try:
                result: T | ErrorOutput | None = input_queue.get()
            except EOFError:
                continue
            if result is None:
                none_counter += 1
                if none_counter == max_nones:
                    # Everything is done, pack up all the results and send them back
                    returnlist: list[T | ErrorOutput] = flatten(
                        list(output_dict.values())
                    )
                    pipe.send(returnlist)
                    logger.debug(
                        "All algorithm returns have been processed, listener is exiting"
                    )
                    break
                continue
            output_dict[result.id].append(result)
            if isinstance(result, ErrorOutput):
                result.write_results()
            if len(output_dict[result.id]) == calls_per_input:
                merged = typer.merge(
                    [
                        x
                        for x in output_dict[result.id]
                        if not isinstance(x, ErrorOutput)
                    ]
                )  # type:ignore Can be fixed with protocol but it gets too complicated and not worth it
                if no_print:
                    continue
                else:
                    for element in merged:
                        writing_started = element.write_results(
                            writing_started, output_file_path
                        )

    def auto_run(
        self,
        inputs: list[ReadyAlgorithmInput],
        no_print: bool,
    ) -> AlgorithmOutput:
        if len(inputs) == 0:
            raise RuntimeError("No input objects provided to auto run function")
        if (
            len(inputs) * inputs[0].call_count > 1 and self.max_processes > 1
        ):  # If either isnt we dont do parallel
            logger.debug(
                "Recognized multiple inputs, creating worker pool and running in parallel"
            )
            return AlgorithmOutput.from_list(
                self._multi_run(inputs, no_print),
                self.result_type,
            )
        else:
            logger.debug("Recognized only one input, running in single process mode")
            return AlgorithmOutput.from_list(
                self._single_run(inputs, no_print), self.result_type
            )

    def _single_run(
        self,
        input_objs: list[ReadyAlgorithmInput],
        no_print: bool,
    ) -> list[T | ErrorOutput]:
        """Function for running a single"""
        returnlist: list[T | ErrorOutput] = []
        initiated = False
        for obj in input_objs:
            for call in obj.call_list:
                try:
                    subprocess_output = subprocess.run(
                        call,
                        text=True,
                        capture_output=True,
                        check=True,
                    )
                    answer = self.result_type(
                        obj.id,
                        subprocess_output.stdout,
                        [subprocess_output.stderr],
                    )
                except subprocess.CalledProcessError as error:
                    answer = ErrorOutput(
                        obj.id,
                        " ".join(call),
                        error.stderr,
                        [str(error.returncode)],
                    )
                    logger.error(answer)
                returnlist.append(answer)
                if no_print:
                    pass
                else:
                    initiated = answer.write_results(initiated, self.output_path)
        return returnlist

    def _multi_run(
        self,
        inputs: list[ReadyAlgorithmInput],
        no_print: bool,
    ) -> list[T | ErrorOutput]:
        # Set Up Everything for a multiprocessed run, first make a multiprocessing manager and fill the worker queue with inputs
        manager = multiprocessing.Manager()
        input_q: multiprocessing.Queue[tuple[list[str], str] | None] = manager.Queue()  # type: ignore Because Queue has type Any
        listener_q = manager.Queue()
        for record in inputs:
            for call in record.call_list:
                input_q.put(CallandID(call, record.id))

        # First we set up  a listener with a connection to the main process
        PipeOut, PipeIn = multiprocessing.Pipe(duplex=False)
        listening = multiprocessing.Process(
            target=self._listener,
            args=(
                listener_q,
                self.output_path,
                PipeIn,
                inputs[0].call_count,
                no_print,
                self.max_processes,
                self.result_type,
            ),
        )
        listening.start()

        # Now we make a pool of workers and
        pool = multiprocessing.Pool(self.max_processes)
        workers: list[multiprocessing.pool.AsyncResult[Any]] = []
        for _ in range(
            self.max_processes
        ):  # Populate the pool with worker functions, each doing nothing but getting items from the input queue and processing them
            work = pool.apply_async(
                SubprocessHandler._worker,
                (input_q, listener_q, self.result_type),
            )
            workers.append(work)
            input_q.put(None)

        # Close the Pool
        pool.close()
        pool.join()
        # Tell the listener that no more objects will be put into the queue with the None sentinel and let him finish working
        listening.join()
        listener_output: list[T | ErrorOutput] = (
            PipeOut.recv()
        )  # Receive the list of outputs from the listener
        return listener_output
