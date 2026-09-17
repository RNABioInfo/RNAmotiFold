import multiprocessing
import multiprocessing.connection
import multiprocessing.pool
import subprocess
from typing import Any, Literal
from pathlib import Path
import src.results.algorithm_output
from src.bgap_rna.input_handler import algorithm_input
from contextlib import redirect_stdout
import sys


class subprocess_handler:
    """Class to manage subprocesses for running RNAmotiFold algorithms."""

    @staticmethod
    def _worker(
        input_queue: "multiprocessing.Queue[algorithm_input]",
        output_queue: "multiprocessing.Queue[src.results.algorithm_output.algorithm_output|src.results.algorithm_output.error]",
        process_type: Literal["mfe", "pfc", "ali"],
    ):
        """Simplest worker function that should work universally, do all pre/post processing outside of this."""
        while not input_queue.empty():
            try:
                input_obj: algorithm_input = input_queue.get_nowait()
            except EOFError:
                break
            subprocess_output = subprocess.run(
                input_obj.runtime_call, text=True, capture_output=True, shell=True
            )
            if subprocess_output.returncode == 0:
                result = src.results.algorithm_output.algorithm_output(
                    input_obj.id, subprocess_output.stdout, [subprocess_output.stderr], process_type
                )
            else:
                result = src.results.algorithm_output.error(input_obj.id, subprocess_output.stderr)
            output_queue.put(result)

    @staticmethod
    def _listener(
        input_queue: "multiprocessing.Queue[src.results.algorithm_output.algorithm_output|src.results.algorithm_output.error|None]",
        output_file: Path | None,
        pipe: multiprocessing.connection.Connection,
    ):
        """Simplest listener funtion that should also work universally, takes the result objects put into its queue by the workers, writes them down or prints them.
        When all workers are done, signaled by the sentinel None in the Queue which comes from the main process, terminates and sends a list of result objects to back.
        """
        return_list: list[
            src.results.algorithm_output.algorithm_output | src.results.algorithm_output.error
        ] = []
        writing_started = False
        while True:
            try:
                result: (
                    src.results.algorithm_output.algorithm_output
                    | src.results.algorithm_output.error
                    | None
                ) = input_queue.get()
            except EOFError:
                continue
            if result is None:
                pipe.send(return_list)
                break
            else:
                if isinstance(result, src.results.algorithm_output.algorithm_output):
                    if isinstance(output_file, Path):
                        with open(output_file, "a+") as write_file:
                            with redirect_stdout(write_file):
                                writing_started = result.write_results(writing_started)
                    else:
                        writing_started = result.write_results(writing_started)
                        sys.stdout.flush()
                return_list.append(result)

    def __init__(
        self,
        max_processes: int,
        inputs: None | list[algorithm_input],
        output_path: Path | None,
        process_type: Literal["mfe", "pfc", "ali"],
    ):
        # Set up the very basics of a new multiprocessing step, How Many processes can we start, what are our inputs and where are we supposed to write the output.
        # If we come from the full RNAmotiFold pipeline the input is pre-processed and output location is confirmed to work -> Alternate contsructor for checking these ?
        self.max_processes = max_processes
        self.inputs = inputs
        self.output_path = output_path
        self.process_type = process_type

    def run(
        self,
    ) -> list[src.results.algorithm_output.algorithm_output | src.results.algorithm_output.error]:
        if self.inputs is None:
            raise ValueError("No inputs provided to subprocess handler")
        # Set Up Everything for a multiprocessed run, first make a multiprocessing manager and fill the worker queue with inputs
        manager = multiprocessing.Manager()
        input_q: multiprocessing.Queue[algorithm_input] = manager.Queue()  # type: ignore Because Queue has type Any
        listener_q = manager.Queue()
        for record in self.inputs:
            input_q.put(record)
        # First we set up  a listener with a connection to the main process
        PipeOut, PipeIn = multiprocessing.Pipe(duplex=False)
        listening = multiprocessing.Process(
            target=self._listener, args=(listener_q, self.output_path, PipeIn)
        )
        listening.start()

        # Now we make a pool of workers and
        pool = multiprocessing.Pool(self.max_processes)
        workers: list[multiprocessing.pool.AsyncResult[Any]] = []
        for _ in range(
            self.max_processes
        ):  # Populate the pool with worker functions, each doing nothing but getting items from the input queue and processing them
            work = pool.apply_async(
                subprocess_handler._worker, (input_q, listener_q, self.process_type)
            )
            workers.append(work)
        # Close the Pool
        pool.close()
        pool.join()
        # Tell the listener that no more objects will be put into the queue with the None sentinel and let him finish working
        listener_q.put(None)
        listening.join()
        listener_output: list[
            src.results.algorithm_output.algorithm_output | src.results.algorithm_output.error
        ] = PipeOut.recv()  # Receive the list of outputs from the listener
        return listener_output
