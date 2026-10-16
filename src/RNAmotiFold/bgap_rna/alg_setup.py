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


class _AlgorithmCompilation:

    def __init__(
        self,
        algorithm: str,
        gapc_path: Path | None = None,
        perl_path: Path | None = None,
        make_path: Path | None = None,
    ):
        self.algorithm: str = algorithm
        self.compiled: bool = False
        self._gapc_path: Path | None = gapc_path
        self._perl_path: Path | None = perl_path
        self._make_path: Path | None = make_path

    def move(self, source_dir: Path, destination_dir: Path):
        logger.debug(
            f"Moving {self.algorithm} from {source_dir} to {destination_dir}"
        )
        source_file: Path = source_dir.joinpath(self.algorithm)
        destination_file: Path = destination_dir.joinpath(
            self.algorithm
        )
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
                options: list[str] = [
                    "-t",
                    "--kbacktrace",
                    "--kbest",
                    "--no-coopt-class",
                ]
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
    def compile_call(self) -> list[str]:
        call: list[str] = []
        call.append(str(self.gapc_path.resolve()))
        call.extend(self.algorithm_options)
        call.append("-i")
        call.append(f"{self.algorithm}")
        call.append("-o")
        call.append(f"{self.algorithm}.cc")
        call.append(self.gap_file)
        return call

    @property
    def add_options_call(self) -> list[str]:
        call: list[str] = []
        call.append(str(self.perl_path.resolve()))
        perlscript_path: str = str(
            Path("Misc") / "Applications" / "addRNAoptions.pl"
        )
        makefile_path = f"{self.algorithm}.mf"
        call.append(perlscript_path)
        call.append(makefile_path)
        call.append("0")
        return call

    @property
    def make_call(self) -> list[str]:
        call: list[str] = []
        call.append(str(self.make_path.resolve()))
        call.append("-f")
        call.append(f"{self.algorithm}.mf")
        return call

    def compile_self(self):
        try:
            transpile_status = subprocess.run(
                self.compile_call,
                check=True,
                cwd=RNAmotiFold.RNALOOPS_PATH,
            )
        except subprocess.CalledProcessError as e:
            raise e
        if transpile_status.returncode == 0:
            try:
                add_options = subprocess.run(
                    self.add_options_call,
                    check=True,
                    cwd=RNAmotiFold.RNALOOPS_PATH,
                )
            except subprocess.CalledProcessError as e:
                raise e
        else:
            raise ValueError(
                f"Transpilation subprocess ended with non-zero status {transpile_status.returncode}"
            )
        if add_options.returncode == 0:
            try:
                making_file = subprocess.run(
                    self.make_call,
                    check=True,
                    cwd=RNAmotiFold.RNALOOPS_PATH,
                )
            except subprocess.CalledProcessError as e:
                raise e
        else:
            raise ValueError(
                f"RNA option adding subprocess ended with non-zero status {transpile_status.returncode}"
            )
        if making_file.returncode == 0:
            filenames = [
                f"{self.algorithm}{x}"
                for x in [
                    ".mf",
                    ".o",
                    ".d",
                    ".cc",
                    ".hh",
                    "_main.cc",
                    "_main.o",
                    "_main.d",
                ]
            ]
            generated_files = [
                RNAmotiFold.RNALOOPS_PATH.joinpath(x).resolve()
                for x in filenames
            ]
            for file in generated_files:
                file.unlink()
            return True
        else:
            raise ValueError(
                f"File maing subprocess ended with non-zero status {making_file.returncode}"
            )

    @property
    def gapc_path(self) -> Path:
        if (
            self._gapc_path is None
        ):  # If nothing was set we check if we can find the dependency
            self._gapc_path = RNAmotiFold.input.dependency_finder.find(
                "gapc"
            )  # This raises Errors if the dependency is not found
        return self._gapc_path

    @property
    def perl_path(self) -> Path:
        if (
            self._perl_path is None
        ):  # If nothing was set we check if we can find the dependency
            self._perl_path = RNAmotiFold.input.dependency_finder.find(
                "perl"
            )  # This raises Errors if the dependency is not found
        return self._perl_path

    @property
    def make_path(self) -> Path:
        if self._make_path is None:
            self._make_path = RNAmotiFold.input.dependency_finder.find(
                "make"
            )
        return self._make_path


def setup_algorithms(
    gapc_path: Path | None,
    perl_path: Path | None,
    make_path: Path | None,
    algorithms: list[str],
    poolboys: int,
) -> bool:
    compilation_list: list[_AlgorithmCompilation] = []
    # Create and list AlgorithmCompilation Objects for each algorithm we want to compile
    for algorithm in algorithms:
        new_obj = _AlgorithmCompilation(
            algorithm, gapc_path, perl_path, make_path
        )
        compilation_list.append(new_obj)

    The_Pool = multiprocessing.Pool(
        processes=poolboys
    )  # Create the pool of workers
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
            comp_obj.move(
                RNAmotiFold.RNALOOPS_PATH, RNAmotiFold.RNAMOTIFOLD_BIN
            )
        return True
    return False


def work_func(comp_obj: _AlgorithmCompilation):
    try:
        comp_obj.compile_self()
        return True
    except subprocess.CalledProcessError as error:
        raise error


def get_dependency(
    configurer: configparser.ConfigParser,
    dependency: str,
) -> Path | None:
    path = configurer.get(
        configurer.default_section, f"{dependency}_path", fallback=None
    )
    if path is None or path == "":
        return RNAmotiFold.input.dependency_finder.find(dependency)
    else:
        return Path(
            configurer.get(
                configurer.default_section, f"{dependency}_path"
            )
        )


def main(
    algorithms: list[str],
    gapc_path: Path | None,
    perl_path: Path | None,
    make_path: Path | None,
    workers: int,
):
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (RNAmotiFold, RNAmoSh, RNAmotiCes)"""
    config = configparser.ConfigParser(allow_no_value=True)
    with open(RNAmotiFold.RNAMOTIFOLD_INSTALL_CONFIG, "r+") as of:
        config.read_file(of, source=str(RNAmotiFold.RNAMOTIFOLD_INSTALL_CONFIG))
    # Ensure we have some GAPC path set somewhere
    if gapc_path is None:
        gapc_path = get_dependency(
            config, "gapc"
        )

    if perl_path is None:
        perl_path = get_dependency(
            config, "perl"
        )

    if make_path is None:
        make_path = get_dependency(
            config, "make"
        )

    logger.info(
        f"Using gapc at {gapc_path}, perl at {perl_path}, and make at {make_path}"
    )
    done = setup_algorithms(gapc_path, perl_path, make_path, algorithms, poolboys=workers)
    if done:
        config.set(config.default_section, "gapc_path", str(gapc_path))
        config.set(config.default_section, "perl_path", str(perl_path))
        config.set(config.default_section, "make_path", str(make_path))
        with open(RNAmotiFold.RNAMOTIFOLD_INSTALL_CONFIG, "w+") as of:
            config.write(of)
        logger.info(
            "Algorithms are all set up, you can now use RNAmotiFold"
        )
    else:
        logger.critical(
            "Something went wrong compiling the RNAmotiFold algorithms, please check outputs"
        )
