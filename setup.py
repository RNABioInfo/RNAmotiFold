import shutil
from pathlib import Path
import subprocess
import argparse
import sys
import src.update_motif_collection

ROOT_DIR = Path(__file__).absolute().parent


def _check_preinstalled_gapc() -> Path:
    """Contains cmd_argument parsing solely for the purpose of checking if an already installed gapc is given"""
    parser = argparse.ArgumentParser(
        prog="SetUp.py",
        description="Set up script for RNAmotiFold. Checks if a modified Bellman's GAP compiler is installed and prepares algorithms.",
        epilog="Does anyone read these anyways?",
    )
    parser.add_argument(
        "preinstalled_gapc_path",
        nargs="?",
        action=preinstalled_check,
        type=Path,
        help="You may input the absolute path to a preinstalled gapcM version. If you don't the script will check if there is already a gapc installed (globally or locally) and if it isn't it will run a CMake Script to set it up locally.",
    )
    args = parser.parse_known_args()
    return args[0].preinstalled_gapc_path


class preinstalled_check(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        if value is None:
            setattr(namespace, self.dest, value)
        else:
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


# Checks if gapc is installed with which or locally
def _detect_gapc() -> str | None:
    """Checks for a gapc installation with which and globs RNAmotiFold folder for any gapc instance (which is presumed to be a modified gapc, if you have a different gapc in here that's on you)"""
    global_gapc = shutil.which("gapc")
    if global_gapc is not None:
        return global_gapc
    else:
        local_gapc = list(ROOT_DIR.glob("**/gapc"))
        try:
            return local_gapc[0]
        except IndexError:
            return None


def setup_algorithms(gapc_path: str):
    RNALOOPS_PATH = _check_submodule("RNALoops")
    RNAMOTIFOLD_BIN = Path.joinpath(ROOT_DIR, "Build", "bin")
    RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
    PERL_PATH = Path.joinpath(RNALOOPS_PATH, "Misc", "Applications", "addRNAoptions.pl")
    compilation_list = []
    for algorithm in ["motmfepretty", "motshapeX", "mothishape_h", "mothishape_b", "mothishape_m"]:
        compilation = f"cd {RNALOOPS_PATH} && {gapc_path} -o {algorithm}.cc -t --kbacktrace --kbest -i {algorithm} RNALoops.gap && perl {PERL_PATH} {algorithm}.mf 0 && make -f {algorithm}.mf && mv {algorithm} {RNAMOTIFOLD_BIN} && rm {algorithm}.o && rm {algorithm}.mf && rm {algorithm}.hh && rm {algorithm}.d && rm {algorithm}.cc && rm {algorithm}_main.o && rm {algorithm}_main.d && rm {algorithm}_main.cc"
        compilation_list.append(compilation)
    for algorithm_subopt_pfc in [
        "motmfepretty_subopt",
        "motshapeX_subopt",
        "mothishape_h_subopt",
        "mothishape_b_subopt",
        "mothishape_m_subopt",
        "motpfc",
        "motshapeX_pfc",
        "mothishape_h_pfc",
        "mothishape_b_pfc",
        "mothishape_m_pfc",
    ]:
        compilation = f"cd {RNALOOPS_PATH} && {gapc_path} -o {algorithm_subopt_pfc}.cc -t -i {algorithm_subopt_pfc} RNALoops.gap && perl {PERL_PATH} {algorithm_subopt_pfc}.mf 0 && make -f {algorithm_subopt_pfc}.mf && mv {algorithm_subopt_pfc} {RNAMOTIFOLD_BIN} && rm {algorithm_subopt_pfc}.o && rm {algorithm_subopt_pfc}.mf && rm {algorithm_subopt_pfc}.hh && rm {algorithm_subopt_pfc}.d && rm {algorithm_subopt_pfc}.cc && rm {algorithm_subopt_pfc}_main.o && rm {algorithm_subopt_pfc}_main.d && rm {algorithm_subopt_pfc}_main.cc"
        compilation_list.append(compilation)
    for comp in compilation_list:
        try:
            subprocess.run(comp, shell=True, check=True)
        except subprocess.CalledProcessError as error:
            raise error
    return True


def _check_submodule(submodule: str) -> Path:
    SUBMOD_DIR = Path.joinpath(ROOT_DIR, "submodules", f"{submodule}")
    if len(list(SUBMOD_DIR.glob("*"))) == 0:
        raise ModuleNotFoundError(
            f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {ROOT_DIR}"
        )
    else:
        return SUBMOD_DIR


def run_cmake():
    print(
        "No preinstalled gapc given or found, building the Bellman's GAP compiler from scratch (This might take a while)"
    )
    BUILD_PATH = Path.joinpath(ROOT_DIR, "Build")
    BUILD_PATH.mkdir(exist_ok=True)
    try:
        build_process = subprocess.run(
            f"cmake ..",
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
            f"cmake --build .",
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


def update_sequences_algorithms():
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (motif, MoSh, MotiCes)"""
    print("Setting up RNAmotiFold, updating RNA 3D Motif sequences...")
    src.update_motif_collection.main()  # fetches latest motif versions
    print("Update finished.\n Checking for commandline arguments for preinstalled gap compiler...")
    preinstalled_gapc_path = _check_preinstalled_gapc()

    if preinstalled_gapc_path is None:
        print(
            "No preinstalled gap compiler set in commandline, checking with which and globbing RNAmotiFold folder..."
        )
        auto_gapc_path = _detect_gapc()
        if auto_gapc_path is None:
            print("No installed gapc found, installing...")
            cmake_generated_gapc_path = run_cmake()
            print("gap compiler installed, installing algorithms...")
            setup_algorithms(cmake_generated_gapc_path)
            print("Algorithms are all set up, you can now use RNAmotiFold")
        else:
            print(f"gap compiler found in {auto_gapc_path}. Using it to set up algorithms...")
            setup_algorithms(auto_gapc_path)
            print("Algorithms are all set up, you can now use RNAmotiFold")
    else:
        print("Preinstalled gap compiler given, using it to install RNAmotiFold...")
        setup_algorithms(preinstalled_gapc_path)
        print("Algorithms are all set up, you can now use RNAmotiFold")


# update function for loading new motif sets into algrithms, DOES NOT LOAD FRESH MOTIF SETS, DOES NOT UPDATE duplicates.json
def update_algorithms():
    """main update function that checks for a gap compiler, installs it if necessary, and (re)compiles all preset algorithms. Does not update motif sequences"""
    src.update_motif_collection.update_hexdumbs()
    preinstalled_gapc_path = _check_preinstalled_gapc()
    if preinstalled_gapc_path is None:
        auto_gapc_path = _detect_gapc()
        if auto_gapc_path is None:
            cmake_generated_gapc_path = run_cmake()
            setup_algorithms(cmake_generated_gapc_path)
        else:
            setup_algorithms(auto_gapc_path)
    else:
        setup_algorithms(preinstalled_gapc_path)


if __name__ == "__main__":
    update_sequences_algorithms()
