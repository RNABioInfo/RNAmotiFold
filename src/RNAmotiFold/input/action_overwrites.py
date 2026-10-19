import argparse
import logging
import subprocess
from collections.abc import Sequence
from os import W_OK, access, cpu_count
from pathlib import Path
from sys import stderr
from typing import Any, Literal

import RNAmotiFold

loggers = logging.getLogger(__name__)

""""Module for different types of check function used during cmd argument parsing and config file parsing. Overwrites of argparse.Action are used to implement
these checks. The check functions are used to check if the given input is valid and if not, raise an error.
"""


class VersionParser(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(namespace, self.dest, VersionParser._parse_version(value))

    @staticmethod
    def _parse_version(value: str | Sequence[Any] | None) -> str:
        if value is not None:
            vers = str(value).replace(".", "_")
            if vers in RNAmotiFold.AVAILABLE_VERSIONS:
                return vers
            else:
                raise ValueError(
                    "Invalid version specified, please select one of the available versions: "
                    + ", ".join(RNAmotiFold.AVAILABLE_VERSIONS)
                )
        else:
            raise ValueError(
                "-v Set but not version specified, please specify a version with -v <version> or --version <version>"
            )


class MotifFileCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace,
            self.dest,
            MotifFileCheck._motif_file_check_function(str(value)),
        )

    @staticmethod
    def _motif_file_check_function(value: str | None) -> Path:
        if value is not None and value != "":
            if Path(value).resolve().is_file():
                return Path(value)
            raise FileNotFoundError("Could not find specified file.")
        raise ValueError("No motif file specified")


class LogCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace,
            self.dest,
            LogCheck._log_check_function(str(value)),
        )

    @staticmethod
    def _log_check_function(value: str) -> str:
        """Checks if the given value is a valid log level, raises ValueError if not"""

        if not isinstance(getattr(logging, value.upper(), None), int):
            raise TypeError(f"Invalid log level: {value}")
        else:
            return value.upper()


class FloatCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace,
            self.dest,
            FloatCheck._float_check_function(str(value)),
        )

    @staticmethod
    def _float_check_function(value: str) -> float:
        """Checks if the given value is a float between 0 and 1, raises ValueError if not"""
        if 1 - float(value) >= 0:
            return float(value)
        else:
            raise ValueError(
                "Invalid Float Value, detected. Please set a value between 0 and 1"
            )


class WorkerCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace, self.dest, WorkerCheck.worker_check_function(value)
        )  # type: ignore

    @staticmethod
    def worker_check_function(
        value: str | Sequence[Any] | None,
    ) -> int:
        """Checks if the given value is a valid number of CPU cores based on the number of available cpus from os.cpu_count, raises ValueError if not"""
        cpus: int | None = cpu_count()
        match value:
            case str():
                if cpus is not None and int(value) > cpus:
                    loggers.info(
                        f"Given worker number exceeds detected cpu count, setting workers to {cpus / 2!s}"
                    )
                    return round(int(cpus) / 2)
                else:
                    return int(value)
            case Sequence() | None:
                loggers.critical(
                    "Couldnt process input number of parallel processes"
                )
                if cpus is not None:
                    loggers.critical(
                        f"Setting number of cpus to {cpus / 2!s}"
                    )
                    return round(int(cpus) / 2)
                else:
                    loggers.critical(
                        "Couldnt access number of CPUs either, playing it save and setting workers to 1"
                    )
                    return 1


class MotifListCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace,
            self.dest,
            MotifListCheck._motif_list_check_function(value),
        )

    @staticmethod
    def _motif_list_check_function(
        value: str | Sequence[Any] | None,
    ) -> str:
        """Checks if the given value is a valid motif list, if input is None or empty, returns an empty string"""
        if value is None:
            return ""
        else:
            return str(value)


class ConfigCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        """Checks if the given value is a valid config file, if input is None or empty, returns an empty string. Raises FileNotFoundError if the file does not exist"""
        if value == "":
            loggers.info("Using default config values")
            setattr(namespace, self.dest, value)
        elif Path(str(value)).resolve().is_file():
            setattr(namespace, self.dest, Path(str(value)))
        else:
            raise FileNotFoundError(
                f"Could not find specified config file {value}"
            )


class OutputFileCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace,
            self.dest,
            OutputFileCheck.output_file_check_function(value),
        )  # type: ignore

    @staticmethod
    def is_path_creatable(pathname: str) -> bool:
        """
        `True` if the current user has sufficient permissions to create the passed
        pathname; `False` otherwise.
        """
        # Parent directory of the passed path. If empty, we substitute the current
        # working directory (CWD) instead.
        dirname = Path(pathname).resolve().parent or Path.cwd()
        return access(dirname, W_OK)

    @staticmethod
    def is_path_exists_or_creatable(pathname: str) -> bool:
        """`True` if the passed pathname is a valid pathname for the current OS _and_
        either currently exists or is hypothetically creatable; `False` otherwise.

        This function is guaranteed to _never_ raise exceptions.
        """
        try:
            # To prevent "os" module calls from raising undesirable exceptions on
            # invalid pathnames, is_pathname_valid() is explicitly called first.
            return (
                Path(pathname).resolve().parent.exists()
                and OutputFileCheck.is_path_creatable(pathname)
            )
        except OSError:
            loggers.error(
                f"Given output file path: {pathname} is neither a file nor a dictionary that the current user can edit, defaulting to outputting to stdout"
            )
            return False

    @staticmethod
    def output_file_check_function(
        value: str | Sequence[Any] | None,
    ) -> str | None:
        """Checks if the given value is a valid output file path, raises FileNotFoundError if not"""
        match value:
            case None:
                return None
            case str():
                if OutputFileCheck.is_path_exists_or_creatable(value):
                    if Path(value).resolve().is_file():
                        stderr.write(
                            f"Given file {value} already exists, results will be appended.\n"
                        )
                    return value
                else:
                    raise FileNotFoundError(
                        "Given path is not a valid path."
                    )
            case Sequence():
                raise TypeError(
                    "Unsupported output file path was invalid of type, trying using a string"
                )


class AlgorithmMatching(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace,
            self.dest,
            AlgorithmMatching.algorithm_matching_function(value),
        )  # type: ignore , ignored cause of the base value typing. Only non protected function

    @staticmethod
    def algorithm_matching_function(
        value: str | Sequence[Any] | None,
    ) -> Literal["RNAmoSh", "RNAmotiCes", "RNAmotiFold", "RNAmotiAlign"]:
        """Checks if the given value is a valid algorithm name, raises ValueError if not"""
        match value:
            case None:
                return "RNAmotiFold"
            case str():
                match value.strip().lower():
                    case "rnamosh":
                        return "RNAmoSh"
                    case "rnamotices":
                        return "RNAmotiCes"
                    case "rnamotifold":
                        return "RNAmotiFold"
                    case "rnamotialign":
                        return "RNAmotiAlign"
                    case _:
                        raise ValueError(
                            f"Invalid algorithm specified: {value}. Valid choices are RNAmoSh, RNAmotiCes, RNAmotiAlign, and RNAmotiFold"
                        )
            case Sequence():
                raise TypeError(
                    "Couldnt process Algorithm input, try inputting a string"
                )


class energyCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace, self.dest, energyCheck._float_check_function(value)
        )

    @staticmethod
    def _float_check_function(value: str | Sequence[Any] | None):
        match value:
            case None:
                raise ValueError(
                    "Can't use energy parameters without setting a value"
                )
            case str():
                if float(value) > 0:
                    return float(value)
            case Sequence():
                raise TypeError("Set energy value is a sequence ?")


class gapcMCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(
            namespace, self.dest, gapcMCheck._gapcm_check_function(value)
        )

    @staticmethod
    def _gapcm_check_function(
        value: str | Sequence[Any] | None,
    ) -> Path | None:
        if isinstance(value, (Path, str)):
            if Path(value).is_file():
                try:
                    version_check = subprocess.run(
                        [f"{value}", "--version"],
                        capture_output=True,
                        check=True,
                    )
                except (
                    subprocess.CalledProcessError,
                    PermissionError,
                ) as error:
                    raise RuntimeError(
                        "Unable to open file, check the above error for more information."
                    ) from error
                else:
                    if "gapcM" in version_check.stdout.decode():
                        return Path(value)
                    else:
                        raise RuntimeError(
                            "The given file is not an instance of the modified Bellman's GAP compiler."
                        )
            else:
                raise FileNotFoundError("The given file does not exist.")
        else:
            raise TypeError("Why is my value a Sequence ?")


class PerlCheck(argparse.Action):
    def __init__(
        self, option_strings: str, dest: str, **kwargs: Any
    ) -> None:
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        setattr(namespace, self.dest, PerlCheckFunction(value))


def PerlCheckFunction(value: str | Sequence[Any] | None) -> Path | None:
    if isinstance(value, str):
        answer = subprocess.run(
            [f"{value}", "-v"], capture_output=True, check=True
        )
        if (
            answer.returncode == 0
            and "This is perl" in answer.stdout.decode()
        ):
            return Path(value).resolve()
        else:
            raise RuntimeError("The given file is not a perl interpreter.")
