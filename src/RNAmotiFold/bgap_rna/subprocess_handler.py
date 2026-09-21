import multiprocessing
import multiprocessing.connection
import multiprocessing.pool
import subprocess
from typing import Any, Literal
from pathlib import Path
import src.RNAmotiFold.results.algorithm_output
from src.RNAmotiFold.bgap_rna.input_handler import algorithm_input
from contextlib import redirect_stdout
import sys
import src.RNAmotiFold.results.mfe

class subprocess_handler:
    """Class to manage subprocesses for running RNAmotiFold algorithms."""

    @staticmethod
    def _worker(
        input_queue: "multiprocessing.Queue[algorithm_input|None]",
        output_queue: "multiprocessing.Queue[src.RNAmotiFold.results.algorithm_output.algorithm_output|src.RNAmotiFold.results.algorithm_output.error]",
        process_type: Literal["mfe", "pfc", "ali"],
    ):
        """Simplest worker function that should work universally, do all pre/post processing outside of this."""
        while True:
            input_obj: algorithm_input|None = input_queue.get()
            if input_obj is None:
                break
            subprocess_output = subprocess.run(
                input_obj.runtime_call, text=True, capture_output=True, shell=True
            )
            if subprocess_output.returncode == 0:
                result = src.RNAmotiFold.results.algorithm_output.algorithm_output(
                    input_obj.id, subprocess_output.stdout, [subprocess_output.stderr], process_type
                )
            else:
                result = src.RNAmotiFold.results.algorithm_output.error(
                    input_obj.id, subprocess_output.stderr
                )
            output_queue.put(result)

    @staticmethod
    def _listener(
        input_queue: "multiprocessing.Queue[src.RNAmotiFold.results.algorithm_output.algorithm_output|src.RNAmotiFold.results.algorithm_output.error|None]",
        output_file: Path | None,
        pipe: multiprocessing.connection.Connection,
        calls_per_input:int,
        merge_mfe:bool,
    ):
        """Simplest listener funtion that should also work universally, takes the result objects put into its queue by the workers, writes them down or prints them.
        When all workers are done, signaled by the sentinel None in the Queue which comes from the main process, terminates and sends a list of result objects to back.
        """
        return_list: list[
            src.RNAmotiFold.results.algorithm_output.algorithm_output
            | src.RNAmotiFold.results.algorithm_output.error
        ] = []
        output_dict:dict[str,list[src.RNAmotiFold.results.algorithm_output.algorithm_output]] = {}
        writing_started = False
        while True:
            try:
                result: (
                    src.RNAmotiFold.results.algorithm_output.algorithm_output
                    | src.RNAmotiFold.results.algorithm_output.error
                    | None
                ) = input_queue.get()
            except EOFError:
                continue
            if result is None:
                pipe.send(return_list)
                break
            else:
                if isinstance(result,src.RNAmotiFold.results.algorithm_output.algorithm_output):
                    output_dict.setdefault(result.id,[]).append(result)
                    if len(output_dict[result.id]) == calls_per_input:
                        match output_dict[result.id][0].process_type:
                            case "mfe":
                                full_output = src.RNAmotiFold.results.algorithm_output.algorithm_output.merge_mfe_outputs(output_dict[result.id])
                                if merge_mfe:
                                    full_output = subprocess_handler.postprocessing_mfe(full_output)
                            case "pfc":
                                full_output = output_dict[result.id]
                                if len(full_output) > 1:
                                    full_output = subprocess_handler.postprocessing_pfc(full_output)
                            case "ali":
                                full_output = output_dict[result.id]
                        if isinstance(output_file, Path):
                            with open(output_file, "a+") as write_file:
                                with redirect_stdout(write_file):
                                    if isinstance(full_output,list):
                                        for element in full_output:
                                            writing_started = element.write_results(writing_started)
                                    else:
                                        writing_started = full_output.write_results(writing_started)
                        else:
                            if isinstance(full_output,list):
                                for element in full_output:
                                     writing_started = element.write_results(writing_started)
                            else:
                                writing_started = full_output.write_results(writing_started)
                            sys.stdout.flush()
                return_list.append(result)

    def __init__(
        self,
        max_processes: int,
        output_path: Path | None,
        process_type: Literal["mfe", "pfc", "ali"],
        calls_per_input:int,
    ):
        # Set up the very basics of a new multiprocessing step, How Many processes can we start, what are our inputs and where are we supposed to write the output.
        # If we come from the full RNAmotiFold pipeline the input is pre-processed and output location is confirmed to work -> Alternate contsructor for checking these ?
        self.max_processes = max_processes
        self.output_path = output_path
        self.process_type = process_type
        self.calls_per_input = calls_per_input

    def run(
        self,inputs:list[algorithm_input],merge_mfe_outputs:bool,
    ) -> list[
        src.RNAmotiFold.results.algorithm_output.algorithm_output
        | src.RNAmotiFold.results.algorithm_output.error
    ]:
        # Set Up Everything for a multiprocessed run, first make a multiprocessing manager and fill the worker queue with inputs
        manager = multiprocessing.Manager()
        input_q: multiprocessing.Queue[algorithm_input|None] = manager.Queue()  # type: ignore Because Queue has type Any
        listener_q = manager.Queue()
        for record in inputs:
            input_q.put(record)
        # First we set up  a listener with a connection to the main process
        PipeOut, PipeIn = multiprocessing.Pipe(duplex=False)
        listening = multiprocessing.Process(
            target=self._listener, args=(listener_q, self.output_path, PipeIn,self.calls_per_input,merge_mfe_outputs)
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
            input_q.put(None)
            
        # Close the Pool
        pool.close()
        pool.join()
        # Tell the listener that no more objects will be put into the queue with the None sentinel and let him finish working
        listener_q.put(None)
        listening.join()
        listener_output: list[
            src.RNAmotiFold.results.algorithm_output.algorithm_output
            | src.RNAmotiFold.results.algorithm_output.error
        ] = PipeOut.recv()  # Receive the list of outputs from the listener
        return listener_output


    #Postprocessing function are not fully implemented yet, update this later FIXME
    @staticmethod
    def postprocessing_pfc(
        merged_output: list[src.RNAmotiFold.results.algorithm_output.algorithm_output],
    ) -> list[src.RNAmotiFold.results.algorithm_output.algorithm_output]:
        returnlist: list[src.RNAmotiFold.results.algorithm_output.algorithm_output] = []
        checklist: list[str] = []
        for output in merged_output:
            if str(output) not in checklist and len(output.results) > 1:
                checklist.append(str(output))
                returnlist.append(output)
        return returnlist

    @staticmethod
    def postprocessing_mfe(merged_output: src.RNAmotiFold.results.algorithm_output.algorithm_output) -> src.RNAmotiFold.results.algorithm_output.algorithm_output:
        """
        Postprocessing function for merging outputs of the seperated motif predictions
        """
        mfe_dict: dict[float, list[src.RNAmotiFold.results.mfe.result_mfe]] = {}
        for res in merged_output.results:
            if (
                isinstance(res, src.RNAmotiFold.results.mfe.result_mfe) and res.classifier != "_"
            ):  # this is a little unnecessary but it gets rid of warnings, the res classifier filter removes the "no motif" structure
                if (
                    res.free_energy not in mfe_dict.keys()
                ):  # -> It makes no sense to have it in the merging process since if it can fit a motif it will be the mfe for that motif anyways
                    mfe_dict[res.free_energy] = [res]
                else:
                    mfe_dict[res.free_energy].append(res)
        for key in mfe_dict.keys():
            if len(mfe_dict[key]) > 1:
                merge_candidates = src.RNAmotiFold.results.mfe.result_mfe.get_compatible_structures(mfe_dict[key])
                for compatible_structures in merge_candidates:
                    new_result = src.RNAmotiFold.results.mfe.result_mfe.merge_structures(
                        [mfe_dict[key][i] for i in compatible_structures]
                    )
                    if new_result is not None:
                        merged_output.results.append(new_result)
            else:
                continue
        merged_output.results.sort(key=lambda x: x.free_energy)  # type: ignore
        return merged_output