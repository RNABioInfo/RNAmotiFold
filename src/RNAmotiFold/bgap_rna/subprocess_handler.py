import multiprocessing
import multiprocessing.connection
import multiprocessing.pool
import subprocess
from typing import Any, Literal
from pathlib import Path
from RNAmotiFold.results.mfe import ResultMFE
from RNAmotiFold.results.algorithm_output import AlgorithmOutput, AlgorithmError
from RNAmotiFold.bgap_rna.input_handler import _AlgorithmInput
from RNAmotiFold.input.parameters import _ScriptParameters
from contextlib import redirect_stdout
import sys
import os
import logging

logger = logging.getLogger(__name__)


class _SubprocessHandler:
    """Class to manage subprocesses for running RNAmotiFold algorithms."""

    @staticmethod
    def _worker(
        input_queue: "multiprocessing.Queue[_AlgorithmInput|None]",
        output_queue: "multiprocessing.Queue[AlgorithmOutput|AlgorithmError]",
        process_type: Literal["mfe", "pfc", "ali"],
    ):
        """Simplest worker function that should work universally, do all pre/post processing outside of this."""
        pid = os.getpid()
        logger.debug(f"Started worker process at {pid}, with process type {process_type}")
        while True:
            input_obj: _AlgorithmInput | None = input_queue.get()
            if input_obj is None:
                logger.debug(f"Worker {pid} exiting, input queue is empty")
                break
            logger.debug(
                f"Worker {pid} started work on {input_obj.id} with call {input_obj.runtime_call}"
            )
            subprocess_output = subprocess.run(
                input_obj.runtime_call, text=True, capture_output=True, shell=True
            )
            if subprocess_output.returncode == 0:
                result = AlgorithmOutput(
                    input_obj.id, subprocess_output.stdout, [subprocess_output.stderr], process_type
                )
                logger.debug(
                    f"Worker {pid} successfully completed call {input_obj.call} on {input_obj.id}"
                )
            else:
                result = AlgorithmError(input_obj.id, subprocess_output.stderr)
                logger.debug(
                    f"Worker {pid} encountered an issue working in {input_obj.id} with call {input_obj.call}"
                )
            output_queue.put(result)

    @staticmethod
    def _listener(
        input_queue: "multiprocessing.Queue[AlgorithmOutput|AlgorithmError|None]",
        output_file: Path | None,
        pipe: multiprocessing.connection.Connection,
        calls_per_input: int,
        merge_mfe: bool,
        no_print: bool,
    ):
        """Simplest listener funtion that should also work universally, takes the result objects put into its queue by the workers, writes them down or prints them.
        When all workers are done, signaled by the sentinel None in the Queue which comes from the main process, terminates and sends a list of result objects to back.
        """
        return_list: list[AlgorithmOutput | AlgorithmError] = []
        output_dict: dict[str, list[AlgorithmOutput]] = {}
        writing_started = False
        while True:
            try:
                result: AlgorithmOutput | AlgorithmError | None = input_queue.get()
            except EOFError:
                continue
            if result is None:
                pipe.send(return_list)
                break
            else:
                if isinstance(result, AlgorithmOutput):
                    return_list.append(result)
                    output_dict.setdefault(result.id, []).append(result)
                    if len(output_dict[result.id]) == calls_per_input:
                        match output_dict[result.id][0].process_type:
                            case "mfe":
                                full_output = AlgorithmOutput._merge_mfe_outputs(
                                    output_dict[result.id]
                                )
                                if merge_mfe:
                                    full_output = _SubprocessHandler.postprocessing_mfe(full_output)
                            case "pfc":
                                full_output = output_dict[result.id]
                                if len(full_output) > 1:
                                    full_output = _SubprocessHandler.postprocessing_pfc(full_output)
                            case "ali":
                                full_output = output_dict[result.id]
                        if isinstance(output_file, Path):
                            with open(output_file, "a+") as write_file:
                                with redirect_stdout(write_file):
                                    if isinstance(full_output, list):
                                        for element in full_output:
                                            writing_started = element.write_results(writing_started)
                                    else:
                                        writing_started = full_output.write_results(writing_started)
                        else:
                            if no_print:
                                continue
                            if isinstance(full_output, list):
                                for element in full_output:
                                    writing_started = element.write_results(writing_started)
                            else:
                                writing_started = full_output.write_results(writing_started)
                            sys.stdout.flush()
                else:
                    logger.critical(
                        f"Error encountered during prediction of {result.id}: {result.error}"
                    )

    def __init__(
        self,
        max_processes: int,
        output_path: Path | None,
        prediction_type: Literal["mfe", "pfc", "ali"],
    ):
        # Set up the very basics of a new multiprocessing step, How Many processes can we start, what are our inputs and where are we supposed to write the output.
        # If we come from the full RNAmotiFold pipeline the input is pre-processed and output location is confirmed to work -> Alternate contsructor for checking these ?
        self.max_processes = max_processes
        self.output_path = output_path
        self.process_type = prediction_type

    @classmethod
    def from_script_parameters(cls, params: _ScriptParameters):
        return cls(params.workers, params.output, params.alg_type())

    def single_run(
        self,
        input_obj: _AlgorithmInput,
        process_type: Literal["mfe", "ali", "pfc"],
        merge_mfe_outputs: bool,
    ) -> AlgorithmOutput | AlgorithmError:
        subprocess_output = subprocess.run(
            input_obj.runtime_call, text=True, capture_output=True, shell=True
        )
        if subprocess_output.returncode == 0:
            result = AlgorithmOutput(
                input_obj.id, subprocess_output.stdout, [subprocess_output.stderr], process_type
            )
            if process_type == "mfe" and merge_mfe_outputs:
                result = self.postprocessing_mfe(result)
            return result
        else:
            return AlgorithmError(input_obj.id, subprocess_output.stderr)

    def run(
        self,
        inputs: list[_AlgorithmInput],
        merge_mfe_outputs: bool,
        calls_per_input: int,
        no_print: bool,
    ) -> list[AlgorithmOutput | AlgorithmError]:
        # Set Up Everything for a multiprocessed run, first make a multiprocessing manager and fill the worker queue with inputs
        manager = multiprocessing.Manager()
        input_q: multiprocessing.Queue[_AlgorithmInput | None] = manager.Queue()  # type: ignore Because Queue has type Any
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
                merge_mfe_outputs,
                no_print,
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
                _SubprocessHandler._worker, (input_q, listener_q, self.process_type)
            )
            workers.append(work)
            input_q.put(None)

        # Close the Pool
        pool.close()
        pool.join()
        # Tell the listener that no more objects will be put into the queue with the None sentinel and let him finish working
        listener_q.put(None)
        listening.join()
        listener_output: list[AlgorithmOutput | AlgorithmError] = (
            PipeOut.recv()
        )  # Receive the list of outputs from the listener
        return listener_output

    @staticmethod
    def postprocessing_pfc(
        merged_output: list[AlgorithmOutput],
    ) -> list[AlgorithmOutput]:
        returnlist: list[AlgorithmOutput] = []
        checklist: list[str] = []
        for output in merged_output:
            if str(output) not in checklist and len(output.results) > 1:
                checklist.append(str(output))
                returnlist.append(output)
        return returnlist

    @staticmethod
    def postprocessing_mfe(merged_output: AlgorithmOutput) -> AlgorithmOutput:
        """
        Postprocessing function for merging outputs of the seperated motif predictions
        """
        mfe_dict: dict[float, list[ResultMFE]] = {}
        for res in merged_output.results:
            if (
                isinstance(res, ResultMFE) and res.classifier != "_"
            ):  # this is a little unnecessary but it gets rid of warnings, the res classifier filter removes the "no motif" structure
                if (
                    res.free_energy not in mfe_dict.keys()
                ):  # -> It makes no sense to have it in the merging process since if it can fit a motif it will be the mfe for that motif anyways
                    mfe_dict[res.free_energy] = [res]
                else:
                    mfe_dict[res.free_energy].append(res)
        for key in mfe_dict.keys():
            if len(mfe_dict[key]) > 1:
                merge_candidates = ResultMFE._get_compatible_structures(mfe_dict[key])
                for compatible_structures in merge_candidates:
                    new_result = ResultMFE._merge_structures(
                        [mfe_dict[key][i] for i in compatible_structures]
                    )
                    if new_result is not None:
                        merged_output.results.append(new_result)
            else:
                continue
        merged_output.results.sort(key=lambda x: x.free_energy)  # type: ignore
        return merged_output
