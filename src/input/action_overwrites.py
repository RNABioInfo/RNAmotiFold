import argparse
from typing import Any, Optional, Sequence, Literal
from sys import stderr
from pathlib import Path
from os import cpu_count, access, W_OK
import logging
from multiprocessing import cpu_count

loggers = logging.getLogger("InputChecks")

""""Module for different types of check function used during cmd argument parsing and config file parsing. Overwrites of argparse.Action are used to implement
these checks. The check functions are used to check if the given input is valid and if not, raise an error.
"""


class MotifFileCheck(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: Optional[str] | Sequence[Any],
        option_string: Optional[str] = None,
    ):
        setattr(namespace, self.dest, MotifFileCheck._motif_file_check_function(str(value)))

    @staticmethod
    def _motif_file_check_function(value: Optional[str]) -> Path:
        if value is not None and value != "":
            if Path(value).resolve().is_file():
                return Path(value)
            raise FileNotFoundError("Could not find specified file.")
        raise ValueError("No motif file specified")


class LogCheck(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: Optional[str] | Sequence[Any],
        option_string: Optional[str] = None,
    ):
        setattr(namespace, self.dest, LogCheck._log_check_function(str(value)))

    @staticmethod
    def _log_check_function(value: str) -> str:
        """Checks if the given value is a valid log level, raises ValueError if not"""        
        if not isinstance(getattr(logging, value.upper(), None), int):
            raise ValueError(f"Invalid log level: {value}")
        else:
            return value.upper()


class FloatCheck(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: Optional[str] | Sequence[Any],
        option_string: Optional[str] = None,
    ):
        setattr(namespace, self.dest, FloatCheck._float_check_function(str(value)))

    @staticmethod
    def _float_check_function(value: str) -> float:
        """Checks if the given value is a float between 0 and 1, raises ValueError if not"""
        if 1.0 - float(value) > 0 or float(value) == 1:
            return float(value)
        else:
            raise ValueError("Invalid Float Value, detected. Please set a value between 0 and 1")


class WorkerCheck(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: Optional[str | Sequence[Any]],
        option_string: Optional[str] = None,
    ):
        setattr(namespace, self.dest, WorkerCheck.worker_check_function(value))  # type: ignore

    @staticmethod
    def worker_check_function(value: Optional[str]) -> Optional[int]:
        """Checks if the given value is a valid number of CPU cores based on the number of available cpus from os.cpu_count, raises ValueError if not"""
        cpus: int = cpu_count()
        if value is not None:
            if int(value) > cpus:
                loggers.info(
                    "Given worker number exceeds detected cpu count, setting workers to cpu_count - 1"
                )
                return int(cpus - 1)
            else:
                return int(value)
        else:
            loggers.info("Could not count cpus, playing it safe and setting CPU_count to 1")
            return 1


class MotifListCheck(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: Optional[str | Sequence[Any]],
        option_string: Optional[str] = None,
    ):
        setattr(namespace, self.dest, MotifListCheck._motif_list_check_function(value))

    @staticmethod
    def _motif_list_check_function(value: Optional[str | Sequence[Any]]) -> str:
        """Checks if the given value is a valid motif list, if input is None or empty, returns an empty string"""
        if value is None:
            return ""
        else:
            return str(value)


class ConfigCheck(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: Optional[str] = None,
    ):
        """Checks if the given value is a valid config file, if input is None or empty, returns an empty string. Raises FileNotFoundError if the file does not exist"""
        if value == "":
            loggers.info(f"Using default config values")
            setattr(namespace, self.dest, value)
        elif Path(str(value)).resolve().is_file():
            setattr(namespace, self.dest, Path(str(value)))
        else:
            raise FileNotFoundError(f"Could not find specified config file {value}")


class OutputFileCheck(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: Optional[str] = None,
    ):
        setattr(namespace, self.dest, OutputFileCheck.output_file_check_function(value))  # type: ignore

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
            return Path(pathname).resolve().parent.exists() and OutputFileCheck.is_path_creatable(
                pathname
            )
        except OSError:
            loggers.error(
                f"Given output file path is neither a file nor a dictionary that the current user can edit, defaulting to outputting to stdout"
            )
            return False

    @staticmethod
    def output_file_check_function(value: Optional[str]):
        """Checks if the given value is a valid output file path, raises FileNotFoundError if not"""
        if value is None:
            return None
        elif OutputFileCheck.is_path_exists_or_creatable(value):
            if Path(value).resolve().is_file():
                stderr.write(f"Given file {value} already exists, results will be appended.\n")
            return value
        else:
            raise FileNotFoundError("Given path is not a valid path.")


class AlgorithmMatching(argparse.Action):
    def __init__(self, option_strings: str, dest: str, **kwargs: Any):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        value: str | Sequence[Any] | None,
        option_string: Optional[str] = None,
    ):
        setattr(namespace, self.dest, AlgorithmMatching.algorithm_matching_function(value))  # type: ignore , ignored cause of the base value typing. Only non protected function

    @staticmethod
    def algorithm_matching_function(
        value: str,
    ) -> Literal["RNAmoSh", "RNAmotiCes", "RNAmotiFold", "RNAmotiAlign"]:
        """Checks if the given value is a valid algorithm name, raises ValueError if not"""
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
