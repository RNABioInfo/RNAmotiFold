import multiprocessing.pool
import shutil
from pathlib import Path
import configparser
import subprocess
import sys
import logging
import multiprocessing
from itertools import product
import RNAmotiFold
import RNAmotiFold.input
import RNAmotiFold.input.cli

ROOT_DIR = Path(__file__).absolute().parents[3]
logger = logging.getLogger(__name__)

class AlgorithmCompilation:

    def __init__(self, algorithm: str, compilescript_call: str):
        self.algorithm: str = algorithm
        self.compilescript_call: str = compilescript_call
        self.compiled: bool = False

    def move(self, source_dir: Path, destination_dir: Path):
        logger.debug(
            f"Moving {self.algorithm} from {source_dir} to {destination_dir}"
        )
        source_file: Path = source_dir.joinpath(self.algorithm)
        destination_file: Path = destination_dir.joinpath(
            self.algorithm
        )
        shutil.move(source_file, destination_file)


def detect_gapc() -> Path:
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
        answer = subprocess.run(
            f"command -v {name}",
            shell=True,
            check=True,
            capture_output=True,
        )
        if (
            answer.returncode == 0
            and f"{name}" in answer.stdout.decode()
        ):
            return Path(answer.stdout.decode()).resolve()
        else:
            raise RuntimeError(
                f"Could not find a {name}, please set path with --{name}_path or install {name} you haven't done so"
            )


def setup_algorithms(
    gapc_path: Path, perl_path: Path, poolboys: int
) -> bool:
    RNALOOPS_PATH = _check_submodule("RNALoops")
    RNAMOTIFOLD_BIN = Path.joinpath(ROOT_DIR, "Build", "bin")
    RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
    COMPILE_SCRIPT = Path.joinpath(
        RNALOOPS_PATH,
        "Misc",
        "Applications",
        "RNAmotiFold",
        "compile.sh",
    )
    compilation_list: list[AlgorithmCompilation] = []
    algorithms = [
        "".join(x)
        for x in list(
            product(
                ["RNAmotiFold", "RNAmoSh", "RNAmotiCes"],
                [
                    "",
                    "Motmicro",
                    "_motmacro_pfc",
                    "_motmacro_subopt",
                    "_subopt",
                    "_pfc",
                ],
            )
        )
    ]
    for algorithm in algorithms:
        if (
            "_" in algorithm
        ):  # There are no motmicro versions of subopt or pfc because of equal structures with different energies in Microstate, see paper Lost in Folding space for details
            options = "-t"
            compilation = f'{COMPILE_SCRIPT} GAPC="{gapc_path}" ALG="{algorithm}" ARGS="{options}" FILE="RNAmotiFold_subopt_pfc.gap" PERL="{perl_path}"'
        else:
            options = "-t --kbacktrace --kbest --no-coopt-class"
            compilation = f'{COMPILE_SCRIPT} GAPC="{gapc_path}" ALG="{algorithm}" ARGS="{options}" FILE="RNAmotiFold.gap" PERL="{perl_path}"'
        compilation_list.append(
            AlgorithmCompilation(algorithm, compilation)
        )

    align = f'{COMPILE_SCRIPT} GAPC="{gapc_path}" ALG="RNAmotiAlign" ARGS="-t --kbacktrace --kbest --no-coopt-class" FILE="RNAmotiAlign.gap" PERL="{perl_path}"'
    compilation_list.append(AlgorithmCompilation("RNAmotiAlign", align))
    if poolboys > len(compilation_list):
        poolboys = len(compilation_list)
    The_Pool = multiprocessing.Pool(processes=poolboys)
    joblist: list[multiprocessing.pool.AsyncResult[bool]] = []
    compilation_success_list: list[bool] = []
    for job in compilation_list:
        obj = The_Pool.apply_async(work_func, (job,))
        joblist.append(obj)
    The_Pool.close()
    The_Pool.join()
    for obj in joblist:
        compilation_success_list.append(obj.successful())
    if all(compilation_success_list):
        logger.info("All compilations successfull, moving binaries...")
        for comp_obj in compilation_list:
            comp_obj.move(RNALOOPS_PATH, RNAMOTIFOLD_BIN)
        return True
    return False


def work_func(comp_obj: AlgorithmCompilation):
    try:
        subprocess.run(
            [comp_obj.compilescript_call], shell=True, check=True
        )
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
        configure_process = subprocess.run(
            f"{cmake_path} ..",
            shell=True,
            check=True,
            stdout=sys.stdout,
            stderr=sys.stdout,
            cwd=BUILD_PATH,
        )
    except subprocess.CalledProcessError as error:
        logger.critical("Error during CMake configuration, exiting...")
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
        logger.critical("Error during CMake building, exiting...")
        raise error
    if (
        not build_process.returncode
        and not configure_process.returncode
    ):
        return Path.joinpath(BUILD_PATH, "gapcM-install", "bin", "gapc")
    raise RuntimeError(
        f"Could not build RNAmotiFold, something went wrong: {build_process.stderr}"
    )


def updates(motif_version: str) -> bool:
    """Does all the updating, fetches perl and gapc paths from defaults or detects them and uses to set up algorithms, returns True if algorithms were updated, False if not"""
    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(
        open(
            file=Path.joinpath(
                ROOT_DIR,
                "src",
                "RNAmotiFold",
                "defaults",
                "defaults.ini",
            )
        )
    )
    update = RNAmotiFold.input.cli.motifs.uninteractive_update(requested_version=motif_version)
    if update:
        if config.get(config.default_section, "perl_path"):
            perl_path = Path(
                config.get(config.default_section, "perl_path")
            )
        else:
            try:
                perl_path = fallback_finder("perl")
            except RuntimeError as error:
                logger.critical(error)
                raise error
        if config.get(config.default_section, "gapc_path"):
            gapc_path = Path(
                config.get(config.default_section, "gapc_path")
            )
        else:
            try:
                gapc_path = detect_gapc()
            except RuntimeError as error:
                logger.critical(error)
                raise error
        if config.get(config.default_section, "setup_workers"):
            poolboys = config.getint(
                config.default_section, "setup_workers"
            )
        else:
            try:
                poolboys = multiprocessing.cpu_count() - 1
            except NotImplementedError as error:
                logger.info(
                    "Could not count cpus, playing it safe and setting CPU_count to 2"
                )
                poolboys = 2
        setup_algorithms(
            gapc_path=gapc_path, perl_path=perl_path, poolboys=poolboys
        )
        return True
    else:
        return False


def main(
    version: str,
    gapc_path: Path | None,
    perl_path: Path | None,
    workers: int,
    cmake_path: Path | None,
):
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (RNAmotiFold, RNAmoSh, RNAmotiCes)"""
    if cmake_path is None:
        cmake_path = fallback_finder("cmake")
    if perl_path is None:
        perl_path = fallback_finder("perl")
    if gapc_path is None:
        try:
            gapc_path = detect_gapc()
        except RuntimeError as error:
            logger.critical(error)
            gapc_path = run_cmake(cmake_path)  # type: ignore

    done: bool = False
    motifs.uninteractive_update(version)  # type: ignore

    done = setup_algorithms(gapc_path, perl_path, workers)
    if done:
        logger.info(
            "Algorithms are all set up, you can now use RNAmotiFold"
        )
    else:
        logger.critical(
            "Something went wrong compiling the RNAmotiFold algorithms, please check outputs"
        )
