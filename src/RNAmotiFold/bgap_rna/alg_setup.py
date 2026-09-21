import multiprocessing.pool
import shutil
from pathlib import Path
import configparser
import subprocess
import argparse
import sys
import logging
import multiprocessing
from itertools import product
import src.RNAmotiFold.input.action_overwrites

ROOT_DIR = Path(__file__).absolute().parents[3]

try:
    import submodules.RNALoops.Misc.Applications.RNAmotiFold.motifs.get_RNA3D_motifs as motifs
except ImportError as e:
    print(
        f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {ROOT_DIR}"
    )
    raise e


logger = logging.getLogger("RNAmotiFold")


def get_cmd_args():
    """Contains cmd_argument parsing solely for the purpose of checking if an already installed gapc is given"""
    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(Path.joinpath(ROOT_DIR, "src","RNAmotiFold", "defaults", "defaults.ini")))
    for option in [
        x for x in config[config.default_section] if config[config.default_section][x] == ""
    ]:
        config.set(config.default_section, option, None)
    parser = argparse.ArgumentParser(
        prog="SetUp.py",
        description="Set up script for RNAmotiFold. Checks if a modified Bellman's GAP compiler is installed and prepares algorithms.",
        epilog="Does anyone read these anyways?",
    )
    parser.add_argument(
        "--cmake_path",
        nargs="?",
        dest="cmake_path",
        action=src.RNAmotiFold.input.action_overwrites.cmake_check,
        default=config.get(config.default_section, "cmake_path"),  # shutil.which("cmake"),
        type=str,
        help=f"Cmake Path for compilation, default can be set at {str(Path.joinpath(ROOT_DIR,"src","RNAmotiFold","defaults","defaults.ini"))}. If no default is set the script will try to find a cmake with which.",
    )
    parser.add_argument(
        "--gapc_path",
        nargs="?",
        action=src.RNAmotiFold.input.action_overwrites.preinstalled_check,
        dest="gapc_path",
        default=config.get(config.default_section, "gapc_path"),  # _detect_gapc(),
        type=str,
        help=f"GAPC Path for compilation, default can be set at {str(Path.joinpath(ROOT_DIR,"src","RNAmotiFold","defaults","defaults.ini"))}.If no default is set the script will try to find a gapc with which and check the RNAmotiFold folder structure for a local installation (it is automatically installed by this script usually).",
    )
    parser.add_argument(
        "--perl_path",
        nargs="?",
        dest="perl_path",
        action=src.RNAmotiFold.input.action_overwrites.perl_check,
        default=config.get(config.default_section, "perl_path"),  # shutil.which("perl"),
        type=str,
        help=f"Perl interpreter path for compilation, default can be set at {str(Path.joinpath(ROOT_DIR,"src","RNAmotiFold","defaults","defaults.ini"))}. If no default is set the script will try to find a perl interpreter with 'which perl' and check /usr/bin/perl.",
    )
    parser.add_argument(
        "-v",
        "--version",
        help=f"Specify which RNA 3D Motif sequence version you want to use. Default is the newest version.",
        dest="version",
        action=src.RNAmotiFold.input.action_overwrites.VersionParser,
        type=str,
        default="current",
    )
    parser.add_argument(
        "-w",
        "-workers",
        type=str,
        dest="workers",
        default=config.get(config.default_section, "setup_workers"),
        help=f"Specify how many parallel processes may be spawned to speed up algorithm compilation. Default can be set at  {str(Path.joinpath(ROOT_DIR,"src","RNAmotiFold","defaults","defaults.ini"))}.",
    )
    args = parser.parse_known_args()[0]

    if args.cmake_path is None:
        cmake_path = fallback_finder("cmake")
        setattr(args, "cmake_path", cmake_path)

    if args.perl_path is None:
        perl_path = fallback_finder("perl")
        setattr(args, "perl_path", perl_path)

    if args.gapc_path is None:
        try:
            gapc_path = _detect_gapc()
        except RuntimeError as error:
            logger.critical(error)
            gapc_path = run_cmake(args.cmake_path)  # type: ignore
        setattr(args, "gapc_path", gapc_path)

    if not args.workers:
        try:
            workers = multiprocessing.cpu_count() - 1
        except NotImplementedError as error:
            logger.critical("Could not count cpus, playing it safe and setting CPU_count to 2")
            workers = 2
        setattr(args, "workers", workers)

    return args

def _detect_gapc() -> Path:
    """Checks for a gapc installation with which and globs RNAmotiFold folder for any gapc instance (which is presumed to be a modified gapc, if you have a different gapc in here that's on you)"""
    global_gapc = shutil.which("gapc")
    if global_gapc is not None:
        return Path(global_gapc)
    else:
        local_gapc = list(ROOT_DIR.glob("**/gapcM-install//bin/gapc"))
        try:
            return local_gapc[0]
        except IndexError:
            raise RuntimeError(
                "Could not find installed gapc, install gapc if necessary or set path to your gapcM executable with --gapc_path or in defaults config"
            )

def fallback_finder(name: str) -> Path:
    whichpath = shutil.which(f"{name}")
    if whichpath is not None:
        return Path(whichpath).resolve()
    else:
        answer = subprocess.run(f"/usr/bin/{name} -v", shell=True, check=True, capture_output=True)
        if answer.returncode == 0 and f"{name}" in answer.stdout.decode():
            return Path(f"/usr/bin/{name}").resolve()
        else:
            raise RuntimeError(
                f"Could not find a {name}, please set path with --{name}_path or install {name} you haven't done so"
            )

def setup_algorithms(gapc_path: Path, perl_path: Path, poolboys: int) -> bool:
    RNALOOPS_PATH = _check_submodule("RNALoops")
    RNAMOTIFOLD_BIN = Path.joinpath(ROOT_DIR, "Build", "bin")
    RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
    COMPILE_SCRIPT = Path.joinpath(
        RNALOOPS_PATH, "Misc", "Applications", "RNAmotiFold", "compile.sh"
    )
    compilation_list: list[str] = []
    algorithms = [
        "".join(x)
        for x in list(
            product(
                ["RNAmotiFold", "RNAmoSh", "RNAmotiCes"],
                ["", "Motmicro", "_motmacro_pfc", "_motmacro_subopt", "_subopt", "_pfc"],
            )
        )
    ]
    for algorithm in algorithms:
        if (
            "_" in algorithm
        ):  # There are no motmicro versions of subopt or pfc because of equal structures with different energies in Microstate, see paper Lost in Folding space for details
            options = "-t"
            compilation = f'{COMPILE_SCRIPT} GAPC="{gapc_path}" ALG="{algorithm}" ARGS="{options}" FILE="RNAmotiFold_subopt_pfc.gap" PERL="{perl_path}" && cd {RNALOOPS_PATH} && mv {algorithm} {RNAMOTIFOLD_BIN}'
        else:
            options = "-t --kbacktrace --kbest --no-coopt-class"
            compilation = f'{COMPILE_SCRIPT} GAPC="{gapc_path}" ALG="{algorithm}" ARGS="{options}" FILE="RNAmotiFold.gap" PERL="{perl_path}" && cd {RNALOOPS_PATH} && mv {algorithm} {RNAMOTIFOLD_BIN}'
        compilation_list.append(compilation)
    align = f'{COMPILE_SCRIPT} GAPC="{gapc_path}" ALG="RNAmotiAlign" ARGS="-t --kbacktrace --kbest --no-coopt-class" FILE="RNAmotiAlign.gap" PERL="{perl_path}" && cd {RNALOOPS_PATH} && mv "RNAmotiAlign" {RNAMOTIFOLD_BIN}'
    The_Pool = multiprocessing.Pool(processes=poolboys)
    compilation_list.append(align)
    joblist: list[multiprocessing.pool.AsyncResult[bool]] = []
    compilation_success_list: list[bool] = []
    for job in compilation_list:
        obj = The_Pool.apply_async(work_func, (job,))
        joblist.append(obj)
    The_Pool.close()
    The_Pool.join()
    for obj in joblist:
        compilation_success_list.append(obj.successful())
    return all(compilation_success_list)

def work_func(call: str):
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

def run_cmake(cmake_path: str | None) -> Path:
    if cmake_path is None:
        raise FileNotFoundError(
            "CMake was not found, please install it or set the path with --cmake_path"
        )
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
        return Path.joinpath(BUILD_PATH, "gapcM-install", "bin", "gapc")
    raise RuntimeError(f"Could not build RNAmotiFold, something went wrong: {build_process.stderr}")

def updates(motif_version: str) -> bool:
    """Does all the updating, fetches perl and gapc paths from defaults or detects them and uses to set up algorithms, returns True if algorithms were updated, False if not"""
    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(file=Path.joinpath(ROOT_DIR, "src","RNAmotiFold", "defaults", "defaults.ini")))
    update = motifs.uninteractive_update(version=motif_version)
    if update:
        if config.get(config.default_section, "perl_path"):
            perl_path = Path(config.get(config.default_section, "perl_path"))
        else:
            try:
                perl_path = fallback_finder("perl")
            except RuntimeError as error:
                logger.critical(error)
                raise error
        if config.get(config.default_section, "gapc_path"):
            gapc_path = Path(config.get(config.default_section, "gapc_path"))
        else:
            try:
                gapc_path = _detect_gapc()
            except RuntimeError as error:
                logger.critical(error)
                raise error
        if config.get(config.default_section, "setup_workers"):
            poolboys = config.getint(config.default_section, "setup_workers")
        else:
            try:
                poolboys = multiprocessing.cpu_count() - 1
            except NotImplementedError as error:
                logger.info("Could not count cpus, playing it safe and setting CPU_count to 2")
                poolboys = 2
        setup_algorithms(gapc_path=gapc_path, perl_path=perl_path, poolboys=poolboys)
        return True
    else:
        return False