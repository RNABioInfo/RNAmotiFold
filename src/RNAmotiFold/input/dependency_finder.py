import shutil
from pathlib import Path
from subprocess import run, CalledProcessError
import logging

logger = logging.getLogger(__name__)


class _dependency_finder:

    def __init__(self, prog: str) -> None:
        self.prog: str = prog

    @staticmethod
    def check_which(prog: str):
        answer = shutil.which(f"{prog}")
        if answer is not None:
            return Path(answer).resolve()
        return None

    @staticmethod
    def check_command(prog: str):
        try:
            answer = run(f"command -v {prog}", shell=True, check=True, capture_output=True)
        except CalledProcessError as e:
            logger.error(f"Error while checking for {prog} with command -v: {e}")
        else:
            if answer.returncode == 0 and f"{prog}" in answer.stdout.decode():
                return Path(answer.stdout.decode()).resolve()
        return None

    def find_dep_path(self):
        dep_path = self.check_which(self.prog)
        if dep_path is None:
            dep_path = self.check_command(self.prog)
        # Add in additional finder steps here if necessary or possible

        if dep_path is not None:
            return dep_path
        raise FileNotFoundError(
            f"Could not find {self.prog} on the system, install if you havent already or set the path to your local installation with --{self.prog}_path"
        )

    @staticmethod
    def check_dependency(prog_path: Path):
        try:
            answer = run([f"{prog_path}", "-v"], capture_output=True, check=True)
        except CalledProcessError as e:
            logger.error(f"Error while checking for {prog_path} with -v: {e}")
            raise e
        else:
            if answer.returncode == 0 and f"{prog_path.name}" in answer.stdout.decode():
                return True
        return False


def _find(dependency: str):
    """Attempts to find a RNAmotiFold dependency with which and command -v /your dependency here/"""
    finder = _dependency_finder(dependency)
    try:
        path = finder.find_dep_path()
        if finder.check_dependency(path):
            return path.resolve()
    except FileNotFoundError as e:
        logger.critical(
            f"Could not find {dependency} on the system, install if you havent already or set it with --{dependency}_path"
        )
        raise e
    except CalledProcessError as e:
        logger.critical(f"Error while checking for {dependency}: {e}")
        raise e
    else:
        logger.critical(
            f"Found dependency {dependency}, but could not confirm it's actually {dependency} with {str(path)} -v, trying to use it anyways"
        )
        return path.resolve()
