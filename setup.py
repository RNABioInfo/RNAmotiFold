import multiprocessing.pool
import shutil
from pathlib import Path
from typing import Optional,Any,Sequence
import subprocess
import argparse
import sys
import logging
import multiprocessing
from itertools import product

ROOT_DIR = Path(__file__).absolute().parent

try:
    import submodules.RNALoops.Misc.Applications.RNAmotiFold.motifs.get_RNA3D_motifs as motifs
except ImportError as e:
    print(
        f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {ROOT_DIR}"
    )


logger = logging.getLogger("RNAmotiFold")


def get_cmd_args():
    """Contains cmd_argument parsing solely for the purpose of checking if an already installed gapc is given"""
    parser = argparse.ArgumentParser(
        prog="SetUp.py",
        description="Set up script for RNAmotiFold. Checks if a modified Bellman's GAP compiler is installed and prepares algorithms.",
        epilog="Does anyone read these anyways?",
    )
    parser.add_argument(
        "--cmake_path",
        nargs="?",
        dest="cmake_path",
        default=shutil.which("cmake"),
        type=str,
        help="If you don't have cmake installed globally or are using a specific CMake version you can set the path here. Default is the cmake version return by 'which cmake'"
    )    
    parser.add_argument(
        "--gapc_path",
        nargs="?",
        action=preinstalled_check,
        dest="preinstalled_gapc_path",
        default=None,
        type=str,
        help="You may input the absolute path to a preinstalled gapcM version. If you don't the script will check if there is already a gapc installed (globally or locally) and if it isn't it will run a CMake Script to set it up locally.",
    )
    parser.add_argument(
        "--perl_path",
        nargs="?",
        dest="perl_path",
        default=shutil.which("perl"),
        type=str,
        help="if you have an alternative perl interpreter you can set the path to it here. Otherwise this script will use the one returned by which perl.",
    )
    parser.add_argument(
        "-v",
        "--version",
        help=f"Specify which RNA 3D Motif sequence version you want to use. Default is the newest version.",
        dest="version",
        type=str,
        default="current",
    )
    parser.add_argument(
        "-w",
        "-workers",
        type=int,
        dest="workers",
        default=5,
        help="Specify how many parallel processes may be spawned to speed up algorithm compilation. Default is 5.",
    )
    args = parser.parse_known_args()
    if args[0].cmake_path is None:
        raise FileNotFoundError("CMake was not found, please install it or set the path with --cmake_path")
    return args[0]


class preinstalled_check(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value: Optional[str|Sequence[Any]], option_string:Optional[str]=None):
        if value is None:
            setattr(namespace, self.dest, value)
        elif isinstance(value,str):
            if Path(value).is_file():
                try:
                    version_check = subprocess.run(
                        [f"{value}", "--version"], capture_output=True, check=True
                    )
                except (subprocess.CalledProcessError, PermissionError) as error:
                    raise RuntimeError(
                        "Unable to open file, check the above error for more information."
                    ) from error
                else:
                    if version_check.stdout.decode()[:4] == "gapc":
                        setattr(namespace, self.dest, Path(value))
                    else:
                        raise RuntimeError(
                            "The given file is not an instance of the modified Bellman's GAP compiler."
                        )
            else:
                raise FileNotFoundError("The given file does not exist.")
        else:
            raise ValueError("Why is my value a Sequence ?")


# Checks if gapc is installed with which or locally
def _detect_gapc() -> Path|None:
    """Checks for a gapc installation with which and globs RNAmotiFold folder for any gapc instance (which is presumed to be a modified gapc, if you have a different gapc in here that's on you)"""
    global_gapc = shutil.which("gapc")
    if global_gapc is not None:
        return Path(global_gapc)
    else:
        local_gapc = list(ROOT_DIR.glob("**/bin/gapc"))
        try:
            return local_gapc[0]
        except IndexError:
            return None

def setup_algorithms(gapc_path: Path, perl_path: Path, poolboys: int) -> bool:
    RNALOOPS_PATH = _check_submodule("RNALoops")
    RNAMOTIFOLD_BIN = Path.joinpath(ROOT_DIR, "Build", "bin")
    RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
    compilation_list:list[str] = []
    algorithms = [
        "".join(x)
        for x in list(product(["RNAmotiFold", "RNAmoSh", "RNAmotiCes"], ["", "_subopt", "_pfc"]))
    ]
    for algorithm in algorithms:
        if "_" in algorithm:
            options = "-t"
            compilation = f'{Path.joinpath(RNALOOPS_PATH,"Misc","Applications","RNAmotiFold","compile.sh")} GAPC="{gapc_path}" ALG="{algorithm}" ARGS="{options}" FILE="RNAmotiFold_subopt_pfc.gap" PERL="{perl_path}" && cd {RNALOOPS_PATH} && mv {algorithm} {RNAMOTIFOLD_BIN}'
        else:
            options = "-t --kbacktrace --kbest --no-coopt-class"
            compilation = f'{Path.joinpath(RNALOOPS_PATH,"Misc","Applications","RNAmotiFold","compile.sh")} GAPC="{gapc_path}" ALG="{algorithm}" ARGS="{options}" FILE="RNAmotiFold.gap" PERL="{perl_path}" && cd {RNALOOPS_PATH} && mv {algorithm} {RNAMOTIFOLD_BIN}'
        compilation_list.append(compilation)
    The_Pool = multiprocessing.Pool(processes=poolboys)
    joblist:list[multiprocessing.pool.AsyncResult[bool]]=[]
    compilation_success_list:list[bool] = []
    for job in compilation_list:
        obj = The_Pool.apply_async(work_func, (job,))
        joblist.append(obj)
    The_Pool.close()
    The_Pool.join()
    for obj in joblist:
        compilation_success_list.append(obj.successful())
    return all(compilation_success_list)

def work_func(call:str):
    try:
        subprocess.run(call, shell=True, check=True)
        return True
    except subprocess.CalledProcessError as error:
        raise error
    
def _check_submodule(submodule: str) -> Path:
    SUBMOD_DIR = Path.joinpath(ROOT_DIR, "submodules", f"{submodule}")
    if len(list(SUBMOD_DIR.glob("*"))) == 0:
        raise ModuleNotFoundError(
            f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {ROOT_DIR}"
        )
    else:
        return SUBMOD_DIR

def run_cmake(cmake_path:str) -> Path:
    BUILD_PATH = Path.joinpath(ROOT_DIR, "Build")
    BUILD_PATH.mkdir(exist_ok=True)
    try:
        build_process = subprocess.run(
            f"{cmake_path} ..",
            shell=True,
            check=True,
            stdout=sys.stdout,
            stderr=sys.stdout,
            cwd=BUILD_PATH,
        )
    except subprocess.CalledProcessError as error:
        print("Error during CMake configuration, exiting...")
        raise error
    try:
        build_process = subprocess.run(
            f"{cmake_path} --build .",
            shell=True,
            check=True,
            stdout=sys.stdout,
            stderr=sys.stdout,
            cwd=BUILD_PATH,
        )
    except subprocess.CalledProcessError as error:
        print("Error during CMake building, exiting...")
        raise error

    if not build_process.returncode:
        return Path.joinpath(BUILD_PATH, "gapc-prefix", "bin", "gapc")
    raise RuntimeError(f"Could not build RNAmotiFold, something went wrong: {build_process.stderr}")

# Does all the updating with tradeoffs between update algorithms and no update
def updates(RNAmotiFold_paramteres:list[str],motif_version: str, workers: int) -> bool:
    update_parser = argparse.ArgumentParser()
    update_parser.add_argument('--perl_path',default=None,dest='perl')
    update_parser.add_argument('--gapc_path',default=None,dest='gapc')
    namespace = update_parser.parse_args(RNAmotiFold_paramteres)
    update = motifs._uninteractive_update(version=motif_version) #type:ignore
    if update:
        if namespace.perl:
            perl_interpreter = namespace.perl
        else:
            perl_interpreter = shutil.which("perl")
            if perl_interpreter is None:
                print(
                    "Could not find a perl interpreter, please input path to your perl interpreter: ",
                    end="",
                )
            perl_interpreter = Path(input())
        if namespace.gapc:
            gapcM_path = namespace.gapc
        else:
            gapcM_path = _detect_gapc()
            if gapcM_path is None:
                print(
                "Could not find gapc, please enter path to your gapcM executeable: ", end=""
                )
                gapcM_path = Path(input())
        setup_algorithms(perl_path=perl_interpreter, gapc_path=gapcM_path, poolboys=workers)
        return True
    else:
        return False


def main():
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (RNAmotiFold, RNAmoSh, RNAmotiCes)"""
    args = get_cmd_args()
    done:bool=False
    if args.preinstalled_gapc_path is None:
        print(
            "No preinstalled gap compiler set in commandline, checking with which and searching RNAmotiFold folder..."
        )
        auto_gapc_path = _detect_gapc()
        if auto_gapc_path is None:
            print("No installed gapc found, installing...")
            cmake_generated_gapc_path = run_cmake(args.cmake_path) #type:ignore
            print("gap compiler installed, installing algorithms...")
            motifs._uninteractive_update(args.version) #type:ignore
            done=setup_algorithms(cmake_generated_gapc_path, args.perl_path, args.workers)
        else:
            print(f"gap compiler found in {auto_gapc_path}. Using it to set up algorithms...")
            motifs._uninteractive_update(args.version) #type:ignore
            done=setup_algorithms(auto_gapc_path, args.perl_path, args.workers)
    else:
        print("Preinstalled gap compiler given, using it to install RNAmotiFold...")
        motifs._uninteractive_update(args.version) #type:ignore
        done=setup_algorithms(args.preinstalled_gapc_path, args.perl_path, args.workers)
    if done:
        print("Algorithms are all set up, you can now use RNAmotiFold")
    else:
        print("Something went wrong compiling the RNAmotiFold algorithms, please check outputs")

if __name__ == "__main__":
    main()
