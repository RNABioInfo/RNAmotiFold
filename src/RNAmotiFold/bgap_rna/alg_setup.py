import multiprocessing.pool
import shutil
from pathlib import Path
import configparser
import subprocess
import logging
import multiprocessing
from itertools import product
import RNAmotiFold
import RNAmotiFold.input.cli
import RNAmotiFold.input.dependency_finder

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
    RNAMOTIFOLD_BIN = Path.joinpath(
        RNAmotiFold._RNAMOTIFOLD_ROOT_DIR, "bin"
    )
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
    SUBMOD_DIR = Path.joinpath(
        RNAmotiFold._RNAMOTIFOLD_ROOT_DIR, f"{submodule}"
    )
    if len(list(SUBMOD_DIR.glob("*"))) == 0:
        raise ModuleNotFoundError(
            f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {RNAmotiFold._RNAMOTIFOLD_ROOT_DIR}"
        )
    else:
        return SUBMOD_DIR


def get_dependency(configurer:configparser.ConfigParser,configpath:Path,dependency:str,user_input_path:Path|str|None) -> Path|None:
    if user_input_path is not None:
        answer = subprocess.run([f"{user_input_path}","-v"],capture_output=True,check=True)
        if answer.returncode == 0 and f"{dependency}" in answer.stdout.decode():
            configurer.set(configurer.default_section,f"{dependency}_path",str(user_input_path))
            with open(configpath,"w+") as of:
                configurer.write(of)
            return Path(user_input_path)
        else:
            logger.error(f"Set path for {dependency}: {str(user_input_path)} could not be identified as a valid instance of {dependency}. Trying prior input...")   
    default_path = configurer.get(configurer.default_section,f"{dependency}_path")
    if default_path:
        logger.error(f"Using {default_path}")
        return Path(default_path)
    else:
        logger.critical(f"Not viable instance of dependency {dependency} was set, neither in {configpath} nor by the user")
        return None

def main(
    version: str,
    gapc_path: Path | None,
    perl_path: Path | None,
    workers: int,
    force_recompile:bool,
    version_update:bool
):
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (RNAmotiFold, RNAmoSh, RNAmotiCes)"""
    config = configparser.ConfigParser(allow_no_value=True)
    configpath = RNAmotiFold._RNAMOTIFOLD_ROOT_DIR / "configs" /"paths.ini"
    with open(configpath,"r+") as of:
        config.read_file(of,source=str(configpath))

    #Check if we need to update anyways because we're on the wrong motif version
    if version_update:
        update = RNAmotiFold.input.cli.motifs.uninteractive_update(version)
    else:
        update = False

    gapc_path = get_dependency(config, configpath, "gapc",gapc_path)
    if gapc_path is None:
        logger.critical("No valid gapc path was set, trying to find it myself")
        gapc_path = RNAmotiFold.input.dependency_finder.find("gapc")
    perl_path = get_dependency(config,configpath,"perl",perl_path)
    if perl_path is None:
        logger.critical("No valid perl path was set, trying to find it myself")
        perl_path = RNAmotiFold.input.dependency_finder.find("perl")
    logger.info(f"Using gapc at {gapc_path} and perl at {perl_path}")
    if update or not force_recompile:
        done = setup_algorithms(gapc_path, perl_path, workers)
        if done:
            logger.info(
                "Algorithms are all set up, you can now use RNAmotiFold"
            )
        else:
            logger.critical(
                "Something went wrong compiling the RNAmotiFold algorithms, please check outputs"
            )
