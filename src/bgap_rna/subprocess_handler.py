import multiprocessing
import multiprocessing.connection
import multiprocessing.pool
import subprocess
from typing import Any
from pathlib import Path
from ..results import base_result
from input_handler import algorithm_input
from contextlib import redirect_stdout
import sys


class SubprocessManager:
    """Class to manage subprocesses for running RNAmotiFold algorithms."""

#    @staticmethod
#    def algorithm_matching(
#        given_algorithm: str,
#    ) -> Literal["RNAmotiFold", "RNAmotiCes", "RNAmoSh", "RNAmotiAlign"]:
#        """MOVE ME TO BGAP RNA SO I DONT HAVE TO PASS ALL THE INFORMATION TO SUBPROCESS HANDLER"""
#        match given_algorithm.lower():
#            case "rnamotifold":
#                return "RNAmotiFold"
#            case "rnamotices":
#                return "RNAmotiCes"
#            case "rnamosh":
#                return "RNAmoSh"
#            case "rnamotialign":
#                return "RNAmotiAlign"
#            case _:
#                raise ValueError(
#                    f"Given Algorithm Name does not match one of the implemented algorithms: {", ".join(["RNAmotiFold","RNAmotiCes","RNAmoSh","RNAmotiAlign"])}."
#                )
#
#    def check_algorithm(self):
#        filepath = (
#            Path(__file__).resolve().parents[1].joinpath("Build", "bin").joinpath(self.algorithm)
#        )
#        if os.path.isfile(filepath) and os.access(filepath, os.X_OK):
#            return True
#        return False

    @staticmethod
    def _worker(
        input_queue: "multiprocessing.Queue[algorithm_input]",
        output_queue: "multiprocessing.Queue[base_result.algorithm_output|base_result.error]"
    ):
        """Simplest worker function that should work universally, do all pre/post processing outside of this."""
        while not input_queue.empty():
            input_obj: algorithm_input = input_queue.get()
            subprocess_output = subprocess.run(
                 input_obj.input_str, text=True, capture_output=True, shell=True
            )
            if subprocess_output.returncode == 0:
                result = base_result.algorithm_output(
                    input_obj.id, subprocess_output.stdout, [subprocess_output.stderr]
                )
            else:
                result = base_result.error(input_obj.id, subprocess_output.stderr)
            output_queue.put(result)

    @staticmethod
    def _listener(
        input_queue: multiprocessing.Queue[base_result.algorithm_output | base_result.error | None],
        output_file: Path | None,
        pipe: multiprocessing.connection.Connection,
    ):
        """Simplest listener funtion that should also work universally, takes the result objects put into its queue by the workers, writes them down or prints them.
        When all workers are done, signaled by the sentinel None in the Queue which comes from the main process, terminates and sends a list of result objects to back.
        """
        return_list: list[base_result.algorithm_output | base_result.error] = []
        writing_started = False
        while True:
            result: base_result.algorithm_output | base_result.error | None = input_queue.get()
            if result is None:
                pipe.send(return_list)
            else:
                if isinstance(result, base_result.algorithm_output):
                    if isinstance(output_file, Path):
                        with open(output_file, "a+") as write_file:
                            with redirect_stdout(write_file):
                                writing_started = result.write_results(writing_started)
                    else:
                        writing_started = result.write_results(writing_started)
                        sys.stdout.flush()
                return_list.append(result)


    def __init__(self, max_processes: int, inputs: list[algorithm_input], output_path: Path | None):
        #Set up the very basics of a new multiprocessing step, How Many processes can we start, what are our inputs and where are we supposed to write the output.
        #If we come from the full RNAmotiFold pipeline the input is pre-processed and output location is confirmed to work -> Alternate contsructor for checking these ?
        self.max_processes = max_processes
        self.inputs = inputs
        self.output_path = output_path


    def run(self) -> list[base_result.algorithm_output | base_result.error]:
        #Set Up Everything for a multiprocessed run, first make a multiprocessing manager and fill the worker queue with inputs
        manager = multiprocessing.Manager()
        input_q: multiprocessing.Queue[algorithm_input] = manager.Queue()  # type: ignore Because Queue has type Any
        listener_q = manager.Queue()
        for record in self.inputs:
            input_q.put(record)
        #First we set up  a listener with a connection to the main process
        PipeOut,PipeIn = multiprocessing.Pipe(duplex=False)
        listening = multiprocessing.Process(target=self._listener,args=(listener_q,self.output_path,PipeIn))
        listening.start()
        
        #Now we make a pool of workers and 
        pool = multiprocessing.Pool(self.max_processes)
        workers:list[multiprocessing.pool.AsyncResult[Any]] = []
        for _ in range(self.max_processes): #Populate the pool with worker functions, each doing nothing but getting items from the input queue and processing them
            work = pool.apply_async(SubprocessManager._worker,(input_q,listener_q))
            workers.append(work)
            
        #This lets the workers finish their all the tasks, no moving past this point until they are all done. Nothing get collected because workers have no return value
        for worker in workers:            
            worker.get()

        #Close the Pool
        pool.close()
        pool.join()
        #Tell the listener that no more objects will be put into the queue with the None sentinel and let him finish working
        listener_q.put(None)
        listening.join()
        listener_output:list[base_result.algorithm_output|base_result.error] = PipeOut.recv() #Receive the list of outputs from the listener 
        return listener_output