import multiprocessing.pool
import shutil
from pathlib import Path
import configparser
import subprocess
import logging
import multiprocessing
import RNAmotiFold
import RNAmotiFold.input.dependency_finder

logger = logging.getLogger(__name__)


class AlgorithmCompilation:

    def __init__(
        self,
        algorithm: str,
        gapc_path: Path | None = None,
        perl_path: Path | None = None,
    ):
        self.algorithm: str = algorithm
        self.compiled: bool = False
        self._gapc_path: Path | None = gapc_path
        self._perl_path: Path | None = perl_path

    def move(self, source_dir: Path, destination_dir: Path):
        logger.debug(f"Moving {self.algorithm} from {source_dir} to {destination_dir}")
        source_file: Path = source_dir.joinpath(self.algorithm)
        destination_file: Path = destination_dir.joinpath(self.algorithm)
        shutil.move(source_file, destination_file)

    @property
    def algorithm_options(self) -> list[str]:
        match self.algorithm:
            case (
                "RNAmotiFold"
                | "RNAmoSh"
                | "RNAmotiCes"
                | "RNAmotiFoldMotmicro"
                | "RNAmoShMotmicro"
                | "RNAmotiCesMotmicro"
                | "RNAmotiAlign"
            ):
                options: list[str] = ["-t", "--kbacktrace", "--kbest", "--no-coopt-class"]
            case (
                "RNAmotiFold_pfc"
                | "RNAmoSh_pfc"
                | "RNAmotiCes_pfc"
                | "RNAmotiFold_motmacro_pfc"
                | "RNAmoSh_motmacro_pfc"
                | "RNAmotiCes_motmacro_pfc"
                | "RNAmotiFold_subopt"
                | "RNAmoSh_subopt"
                | "RNAmotiCes_subopt"
                | "RNAmotiFold_motmacro_pfc"
                | "RNAmotiFold_motmacro_subopt"
                | "RNAmoSh_motmacro_subopt"
                | "RNAmotiCes_motmacro_subopt"
            ):
                options: list[str] = ["-t"]
            case _:
                raise ValueError(
                    f"Unrecognized algorithm {self.algorithm}, could not set options for compilation"
                )
        return options

    @property
    def gap_file(self) -> str:
        if "pfc" in self.algorithm or "subopt" in self.algorithm:
            return "RNAmotiFold_subopt_pfc.gap"
        elif self.algorithm == "RNAmotiAlign":
            return "RNAmotiAlign.gap"
        else:
            return "RNAmotiFold.gap"

    @property
    def compile_call(self) -> str:

        COMPILE_SCRIPT = Path.joinpath(
            RNAmotiFold.RNALOOPS_PATH,
            "Misc",
            "Applications",
            "RNAmotiFold",
            "compile.sh",
        ).resolve()
        return f'{COMPILE_SCRIPT} GAPC="{self.gapc_path}" ALG="{self.algorithm}" ARGS="{" ".join(self.algorithm_options)}" FILE="{self.gap_file}" PERL="{self.perl_path}"'

    @property
    def gapc_path(self) -> Path:
        if self._gapc_path is None:  # If nothing was set we check if we can find the dependency
            self._gapc_path = RNAmotiFold.input.dependency_finder.find(
                "gapc"
            )  # This raises Errors if the dependency is not found
        return self._gapc_path

    @property
    def perl_path(self) -> Path:
        if self._perl_path is None:  # If nothing was set we check if we can find the dependency
            self._perl_path = RNAmotiFold.input.dependency_finder.find(
                "perl"
            )  # This raises Errors if the dependency is not found
        return self._perl_path


def setup_algorithms(
    gapc_path: Path | None, perl_path: Path | None, algorithms: list[str], poolboys: int
) -> bool:
    compilation_list: list[AlgorithmCompilation] = []
    # Create and list AlgorithmCompilation Objects for each algorithm we want to compile
    for algorithm in algorithms:
        new_obj = AlgorithmCompilation(algorithm, gapc_path, perl_path)
        compilation_list.append(new_obj)

    if poolboys > len(
        compilation_list
    ):  # If we need less workers than available we set the number of workers to the number of algorithms to compile
        poolboys = len(compilation_list)

    The_Pool = multiprocessing.Pool(processes=poolboys)  # Create the pool of workers
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
            comp_obj.move(RNAmotiFold.RNALOOPS_PATH, RNAmotiFold.RNAMOTIFOLD_BIN)
        return True
    return False


def work_func(comp_obj: AlgorithmCompilation):
    try:
        subprocess.run([comp_obj.compile_call], shell=True, check=True)
        return True
    except subprocess.CalledProcessError as error:
        raise error


def get_dependency(
    configurer: configparser.ConfigParser, configpath: Path, dependency: str
) -> Path | None:
    path = configurer.get(configurer.default_section, f"{dependency}_path", fallback=None)
    if path is None or path == "":
        return RNAmotiFold.input.dependency_finder.find(dependency)
    else:
        return Path(configurer.get(configurer.default_section, f"{dependency}_path"))


def main(
    algorithms: list[str],
    gapc_path: Path | None,
    perl_path: Path | None,
    workers: int,
):
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (RNAmotiFold, RNAmoSh, RNAmotiCes)"""
    config = configparser.ConfigParser(allow_no_value=True)
    with open(RNAmotiFold.RNAMOTIFOLD_PATHS_CONFIG, "r+") as of:
        config.read_file(of, source=str(RNAmotiFold.RNAMOTIFOLD_PATHS_CONFIG))
    # Ensure we have some GAPC path set somewhere
    if gapc_path is None:
        gapc_path = get_dependency(config, RNAmotiFold.RNAMOTIFOLD_PATHS_CONFIG, "gapc")

    if perl_path is None:
        perl_path = get_dependency(config, RNAmotiFold.RNAMOTIFOLD_PATHS_CONFIG, "perl")

    logger.info(f"Using gapc at {gapc_path} and perl at {perl_path}")
    done = setup_algorithms(gapc_path, perl_path, algorithms, workers)
    if done:
        config.set(config.default_section, "gapc_path", str(gapc_path))
        config.set(config.default_section, "perl_path", str(perl_path))
        with open(RNAmotiFold.RNAMOTIFOLD_PATHS_CONFIG, "w+") as of:
            config.write(of)
        logger.info("Algorithms are all set up, you can now use RNAmotiFold")
    else:
        logger.critical(
            "Something went wrong compiling the RNAmotiFold algorithms, please check outputs"
        )
