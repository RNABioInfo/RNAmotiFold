import argparse
import subprocess
import configparser
from typing import Optional
import sys
import logging
import motif_collection as mc
import args
from pathlib import Path
from dataclasses import dataclass
import process


@dataclass
class script_parameters:
    loglevel: str
    workers: int
    separator: str
    force_update: bool
    no_update: bool
    remove_bool: bool
    pfc_filtering: bool
    custom_algorithm_comp: Optional[str] = None

    @classmethod
    def from_argparse(cls, cmd_args: argparse.Namespace):
        return cls(
            cmd_args.loglevel,
            cmd_args.workers,
            cmd_args.separator,
            cmd_args.force_update,
            cmd_args.no_update,
            cmd_args.remove_bool,
            cmd_args.pfc_filtering_bool,
        )

    @classmethod
    def from_config(cls, config_path: str | Path | configparser.ConfigParser):
        if isinstance(config_path, str) or isinstance(config_path, Path):
            config = configparser.ConfigParser()
            config.read_file(open(config_path))
            return cls(
                config["RNALOOPS"]["loglevel"],
                config.getint("RNALOOPS", "workers"),
                config["RNALOOPS"]["separator"],
                config.getboolean("RNALOOPS", "force_update"),
                config.getboolean("RNALOOPS", "no_update"),
                config.getboolean("RNALOOPS", "remove_bool"),
                config.getboolean("RNALOOPS", "pfc_filtering_bool"),
            )
        elif isinstance(config_path, configparser.ConfigParser):
            return cls(
                config_path["RNALOOPS"]["loglevel"],
                config_path.getint("RNALOOPS", "workers"),
                config_path["RNALOOPS"]["separator"],
                config_path.getboolean("RNALOOPS", "force_update"),
                config_path.getboolean("RNALOOPS", "no_update"),
                config_path.getboolean("RNALOOPS", "remove_bool"),
                config_path.getboolean("RNALOOPS", "pfc_filtering_bool"),
            )

    @staticmethod
    def get_conf_path() -> Path:
        return Path(__file__).resolve().parent.joinpath("data", "config.ini")

    @staticmethod
    def get_RNALoops_path() -> Path:
        return Path(__file__).resolve().parents[1]

    @staticmethod
    def setup_RNALoops_alorithm(proc_obj: process.bgap_rna) -> Path:
        if script_parameters.get_RNALoops_path().joinpath(proc_obj.algorithm) not in list(
            script_parameters.get_RNALoops_path().glob("*")
        ):
            if proc_obj.pfc:
                compilation = f"cd {script_parameters.get_RNALoops_path()} && gapc -o {proc_obj.algorithm}.cc -t -i {proc_obj.algorithm} RNALoops.gap && perl Misc/Applications/addRNAoptions.pl {proc_obj.algorithm}.mf 0 && make -f {proc_obj.algorithm}.mf"
            elif proc_obj.subopt:
                compilation = f"cd {script_parameters.get_RNALoops_path()} && gapc -o {proc_obj.algorithm}.cc -t --kbacktrace -i {proc_obj.algorithm} RNALoops.gap && perl Misc/Applications/addRNAoptions.pl {proc_obj.algorithm}.mf 0 && make -f {proc_obj.algorithm}.mf"
            else:
                compilation = f"cd {script_parameters.get_RNALoops_path()} && gapc -o {proc_obj.algorithm}.cc -t --no-coopt-class --kbacktrace --kbest -i {proc_obj.algorithm} RNALoops.gap && perl Misc/Applications/addRNAoptions.pl {proc_obj.algorithm}.mf 0 && make -f {proc_obj.algorithm}.mf"
            compilation_return = subprocess.run(
                compilation, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            if compilation_return.returncode:
                raise RuntimeError(compilation_return.stderr)
            else:
                subprocess.run(
                    "cd {path} && rm {alg}.o && rm {alg}.mf && rm {alg}.hh && rm {alg}.d && rm {alg}.cc && rm {alg}_main.o && rm {alg}_main.d && rm {alg}_main.cc && rm string.d && rm string.o".format(
                        path=script_parameters.get_RNALoops_path(), alg=proc_obj.algorithm
                    ),
                    shell=True,
                    capture_output=False,
                )
                path = script_parameters.get_RNALoops_path().joinpath(proc_obj.algorithm)
            return path
        else:
            return script_parameters.get_RNALoops_path().joinpath(proc_obj.algorithm)


# Non Process class function that need to be unbound to be pickle'able. See: https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map, guess it kinda is possible it really isnt all that necessary though.


def make_new_logger(
    lvl: str,
    name: str,
    form: str = "",
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(lvl.upper())
    handler = logging.StreamHandler(sys.stderr)
    if form:
        formatter = logging.Formatter(fmt=form)
    else:
        formatter = logging.Formatter(fmt="%(asctime)s:%(levelname)s: %(message)s")
    handler.setFormatter(formatter)
    logger.propagate = (
        False  # Permanent propagate False since I'm not utilizing any subloggers and sublcasses.
    )
    logger.addHandler(handler)
    return logger


def update_checks(
    proc_obj: process.bgap_rna,
    installed_motif_version: str,
    update_force: bool,
    seq_remove_bool: bool,
    logger: logging.Logger,
) -> str | None:
    logger.info("Checking motif version")
    update = proc_obj.check_motif_versions(installed_motif_version)
    if update:
        logger.info("RNA 3D Motif Atlas motif set is outdated, updating...")
    else:
        logger.info("RNA 3D Motif Atlas motif set is up to date.")
        if update_force:
            logger.info("Force update enabled, updating RNA 3D Motif Atlas motif set.")
    if update or update_force:
        mc.update(logger, seq_remove_bool, script_parameters.get_RNALoops_path())
        logger.info(
            f"Update successful, updated from {installed_motif_version} to {process.bgap_rna.get_current_motif_version()}. Recompiling algorithm."
        )
        return True
    return False


if __name__ == "__main__":
    commandline_arguments = args.get_cmdarguments()[0]
    config = args.get_config(script_parameters.get_conf_path())
    if commandline_arguments.config:
        proc = process.bgap_rna.from_config(script_parameters.get_conf_path())
        params = script_parameters.from_config(config)
    else:
        proc = process.bgap_rna.from_argparse(commandline_arguments)
        params = script_parameters.from_argparse(commandline_arguments)
    logger = make_new_logger(params.loglevel, __name__)
    if not params.no_update:
        update = update_checks(
            proc, config["VERSIONS"]["motifs"], params.force_update, params.remove_bool, logger
        )
        if update:
            config["VERSIONS"]["motifs"] = process.bgap_rna.get_current_motif_version()
            with open(script_parameters.get_conf_path(), "w") as file:
                args.get_config(script_parameters.get_conf_path()).write(file)
    else:
        logger.info("Updating disabled, continuing without updating motifs")
    script_parameters.setup_RNALoops_alorithm(proc)
    proc.auto_run(
        writing=True,
        pool_workers=params.workers,
        result_algorithm_csv_separator=params.separator,
    )
