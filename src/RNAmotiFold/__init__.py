import importlib.util as _importlib
from pathlib import Path

RNAMOTIFOLD_ROOT_DIR = Path(__file__).resolve().absolute().parent


def _check_submodule(submodule: str) -> Path:
    SUBMOD_DIR = Path.joinpath(RNAMOTIFOLD_ROOT_DIR, f"{submodule}")
    if len(list(SUBMOD_DIR.glob("*"))) == 0:
        raise ModuleNotFoundError(
            f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {RNAMOTIFOLD_ROOT_DIR}"
        )
    else:
        return SUBMOD_DIR


try:
    _script_dir: Path = (
        RNAMOTIFOLD_ROOT_DIR
        / "RNALoops"
        / "Misc"
        / "Applications"
        / "RNAmotiFold"
        / "motifs"
        / "get_rna3d_motifs.py"
    )
    _spec = _importlib.spec_from_file_location("uniteractive_update", _script_dir)
    if _spec is None or _spec.loader is None:
        raise ImportError("Submodule RNALoops was not correctly cloned.")
    _motifs = _importlib.module_from_spec(_spec)
    _spec.loader.exec_module(_motifs)
except ImportError as e:
    print(e)


_RNAMOTIFOLD_CONFIG_DIR: Path = Path.joinpath(RNAMOTIFOLD_ROOT_DIR, "configs")
RNAMOTIFOLD_DEFAULTS_CONFIG = Path.joinpath(_RNAMOTIFOLD_CONFIG_DIR, "defaults.ini")
RNAMOTIFOLD_INSTALL_CONFIG: Path = Path.joinpath(_RNAMOTIFOLD_CONFIG_DIR, "install.ini")
RNAMOTIFOLD_BIN: Path = Path.joinpath(RNAMOTIFOLD_ROOT_DIR, "bin")
RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
RNALOOPS_PATH: Path = _check_submodule("RNALoops")
_RNAMOTIFOLD_MOTIFS_PATH: Path = Path.joinpath(
    RNALOOPS_PATH,
    "Misc",
    "Applications",
    "RNAmotiFold",
    "motifs",
    "versions",
    "combined",
)
EMPTY_CSV_PATH = (
    RNAMOTIFOLD_ROOT_DIR
    / "RNALoops"
    / "Misc"
    / "Applications"
    / "RNAmotiFold"
    / "motifs"
    / "versions"
    / "empty.csv"
)

AVAILABLE_VERSIONS: list[str] = [
    x.name for x in _RNAMOTIFOLD_MOTIFS_PATH.iterdir() if x.is_dir()
]

from RNAmotiFold.api import (
    rnamosh,
    rnamosh_pfc,
    rnamosh_subopt,
    rnamotialign,
    rnamotices,
    rnamotices_pfc,
    rnamotices_subopt,
    rnamotifold,
    rnamotifold_pfc,
    rnamotifold_subopt,
)

__all__ = (
    "rnamosh",
    "rnamosh_pfc",
    "rnamosh_subopt",
    "rnamotialign",
    "rnamotices",
    "rnamotices_pfc",
    "rnamotices_subopt",
    "rnamotifold",
    "rnamotifold_pfc",
    "rnamotifold_subopt",
)
