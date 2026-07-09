import multiprocessing.pool
import shutil
from pathlib import Path
from typing import Optional,Any,Sequence
import configparser
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
    print(f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {ROOT_DIR}")
    raise e        


logger = logging.getLogger("RNAmotiFold")


def get_cmd_args():
    """Contains cmd_argument parsing solely for the purpose of checking if an already installed gapc is given"""
    config= configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(Path.joinpath(ROOT_DIR,"src","data","defaults.ini")))
    for option in [x for x in config["SETUP"] if config["SETUP"][x] == ""]:
        config.set("SETUP", option, None)
    parser = argparse.ArgumentParser(
        prog="SetUp.py",
        description="Set up script for RNAmotiFold. Checks if a modified Bellman's GAP compiler is installed and prepares algorithms.",
        epilog="Does anyone read these anyways?",
    )
    parser.add_argument(
        "--cmake_path",
        nargs="?",
        dest="cmake_path",
        default=config.get("SETUP","cmake_path"), #shutil.which("cmake"),
        type=str,
        help=f"Cmake Path for compilation, default can be set at {str(Path.joinpath(ROOT_DIR,'src','data','defaults.ini'))} under 'SETUP'. If nod default is set the script will try to find a cmake with which."
    )    
    parser.add_argument(
        "--gapc_path",
        nargs="?",
        action=preinstalled_check,
        dest="preinstalled_gapc_path",
        default=config.get("SETUP","gapc_path"),#_detect_gapc(),
        type=str,
        help=f"GAPC Path for compilation, default can be set at {str(Path.joinpath(ROOT_DIR,"src","data","defaults.ini"))} under 'SETUP'.If no default is set the script will try to find a gapc with which and check the RNAmotiFold folder structure for a local installation (it is automatically installed by this script usually).",
    )
    parser.add_argument(
        "--perl_path",
        nargs="?",
        dest="perl_path",
        action=perl_check,
        default=config.get("SETUP","perl_path"),#shutil.which("perl"),
        type=str,
        help=f"Perl interpreter path for compilation, default can be set at {str(Path.joinpath(ROOT_DIR,"src","data","defaults.ini"))} under 'SETUP'. If no default is set the script will try to find a perl interpreter with 'which perl' and check /usr/bin/perl.",
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
        default=config.get("SETUP","setup_workers"),
        help=f"Specify how many parallel processes may be spawned to speed up algorithm compilation. Default is {config.get("SETUP","workers"),}.",
    )
    args = parser.parse_known_args()
    return args[0]

class cmake_check(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value: Optional[str|Sequence[Any]], option_string:Optional[str]=None):
        if value is None:
            detected_cmake_path = shutil.which("cmake")
            if detected_cmake_path is None:
                raise FileNotFoundError("CMake was not found, please install it or set the path with --cmake_path")
            setattr(namespace, self.dest, detected_cmake_path)
        else:
            if Path(str(value)).is_file():
                try:
                    version_check = subprocess.run(
                        [f"{value}", "--version"], capture_output=True, check=True
                    )
                except (subprocess.CalledProcessError, PermissionError) as error:
                    raise RuntimeError(
                        "Unable to open file, check the above error for more information."
                    ) from error
                else:
                    if "cmake version" in version_check.stdout.decode().lower():
                        setattr(namespace, self.dest, value)
                    else:
                        raise RuntimeError(
                            "The given file is not an instance of CMake."
                        )

class preinstalled_check(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value: Optional[str|Sequence[Any]], option_string:Optional[str]=None):
        if value is None:
            try:
                value = str(_detect_gapc())
            except RuntimeError as error:
                logger.critical("No gapc detected, trying to install with cmake")
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
                    if "gapc" in version_check.stdout.decode()[:4]:
                        setattr(namespace, self.dest, Path(value))
                    else:
                        raise RuntimeError(
                            "The given file is not an instance of the modified Bellman's GAP compiler."
                        )
            else:
                raise FileNotFoundError("The given file does not exist.")
        else:
            raise ValueError("Why is my value a Sequence ?")

class perl_check(argparse.Action):
    def __init__(self, option_strings:str, dest:str, **kwargs:Any):
        super().__init__(option_strings, dest, **kwargs)
    def __call__(self, parser:argparse.ArgumentParser, namespace:argparse.Namespace, value: Optional[str|Sequence[Any]], option_string:Optional[str]=None):
        setattr(namespace,self.dest,PerlCheckFunction(value))

def PerlCheckFunction(value:Optional[str|Sequence[Any]]) -> Path|None:
    if value is None:
        try:
            return fallback_perl_finder()
        except RuntimeError as error:
            logger.critical(error)
            raise error
    elif isinstance(value,str):
        answer = subprocess.run([f"{value}", "-v"],capture_output=True,check=True)
        if answer.returncode == 0 and "This is perl" in answer.stdout.decode():
            return Path(value).resolve()
        else:
            raise RuntimeError("The given file is not a perl interpreter.")

def _detect_gapc() -> Path:
    """Checks for a gapc installation with which and globs RNAmotiFold folder for any gapc instance (which is presumed to be a modified gapc, if you have a different gapc in here that's on you)"""
    global_gapc = shutil.which("gapc")
    if global_gapc is not None:
        return Path(global_gapc)
    else:
        local_gapc = list(ROOT_DIR.glob("**/bin/gapc"))
        try:
            return local_gapc[0]
        except IndexError:
            raise RuntimeError("Could not find insalled gapc, install gapc if necessary or set path to your gapcM executable with --gapc_path")

def fallback_perl_finder() -> Path:
    whichpath = shutil.which("perl")
    if whichpath is not None:
        return Path(whichpath).resolve()
    else:
        answer = subprocess.run("/usr/bin/perl -v",shell=True,check=True,capture_output=True)
        if answer.returncode == 0 and "This is perl" in answer.stdout.decode():
            return Path("/usr/bin/perl").resolve()
        else:
            raise RuntimeError("Could not find a perl interpreter, please set path to your perl interpreter with --perl_path or install perl you haven't done so")

def setup_algorithms(gapc_path: Path, perl_path: Path, poolboys: int) -> bool:
    RNALOOPS_PATH = _check_submodule("RNALoops")
    RNAMOTIFOLD_BIN = Path.joinpath(ROOT_DIR, "Build", "bin")
    RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
    compilation_list:list[str] = []
    algorithms = [
        "".join(x)
        for x in list(product(["RNAmotiFold", "RNAmoSh", "RNAmotiCes"], ["","Motmicro","_motmacro_pfc","_motmacro_subopt", "_subopt", "_pfc"]))
    ]
    for algorithm in algorithms:
        if "_" in algorithm: #There are no motmicro versions of subopt or pfc because of equal structures with different energies in Microstate, see paper Lost in Folding space for details
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

def run_cmake(cmake_path:Optional[str]) -> Path:
    if cmake_path is None:
        raise FileNotFoundError("CMake was not found, please install it or set the path with --cmake_path")
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

def updates(motif_version: str) -> bool:
    '''Does all the updating, fetches perl and gapc paths from defaults or detects them and uses to set up algorithms, returns True if algorithms were updated, False if not'''
    config= configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(file=Path.joinpath(ROOT_DIR,"src","data","defaults.ini")))

    update = motifs._uninteractive_update(version=motif_version) #type:ignore
    if update:
        if config.get("SETUP","perl_path"):
            perl_path = Path(config.get("SETUP","perl_path"))
        else:
            try:
                perl_path = fallback_perl_finder()
            except RuntimeError as error:
                logger.critical(error)
                raise error
        if config.get("SETUP","gapc_path"):
            gapc_path = Path(config.get("SETUP","gapc_path"))
        else:
            try:
                gapc_path = _detect_gapc()
            except RuntimeError as error:
                logger.critical(error)
                raise error
        if config.get("SETUP","setup_workers"):
            poolboys = config.getint("SETUP","setup_workers")
        else:
            try:
                poolboys = multiprocessing.cpu_count() - 1
            except NotImplementedError as error:
                logger.critical("Could not count cpus, playing it safe and setting CPU_count to 2")
                poolboys = 2
        setup_algorithms(gapc_path=gapc_path, perl_path=perl_path, poolboys=poolboys)
        return True
    else:
        return False

def main():
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (RNAmotiFold, RNAmoSh, RNAmotiCes)"""
    args = get_cmd_args()
    done:bool=False
    if args.preinstalled_gapc_path is None:
        logger.critical("No gapc given or detected, trying to install with cmake")
        gapc_path = run_cmake(args.cmake_path) #type:ignore
        logger.critical("gap compiler installed, setting up algorithms")
        print("gap compiler installed, setting up algorithms")
    else:
        logger.critical(f"Preinstalled gap compiler given as {args.preinstalled_gapc_path}")
        gapc_path = args.preinstalled_gapc_path   
    motifs._uninteractive_update(args.version) #type:ignore
    done=setup_algorithms(gapc_path, args.perl_path, args.workers)
    if done:
        print("Algorithms are all set up, you can now use RNAmotiFold")
    else:
        print("Something went wrong compiling the RNAmotiFold algorithms, please check outputs")

if __name__ == "__main__":
    main()
