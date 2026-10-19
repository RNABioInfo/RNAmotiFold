import importlib.util as _importlib
from pathlib import Path as _Path

RNAMOTIFOLD_ROOT_DIR = _Path(__file__).resolve().absolute().parent


def _check_submodule(submodule: str) -> _Path:
    SUBMOD_DIR = _Path.joinpath(RNAMOTIFOLD_ROOT_DIR, f"{submodule}")
    if len(list(SUBMOD_DIR.glob("*"))) == 0:
        raise ModuleNotFoundError(
            f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {RNAMOTIFOLD_ROOT_DIR}"
        )
    else:
        return SUBMOD_DIR


try:
    _script_dir: _Path = (
        RNAMOTIFOLD_ROOT_DIR
        / "RNALoops"
        / "Misc"
        / "Applications"
        / "RNAmotiFold"
        / "motifs"
        / "get_rna3D_motifs.py"
    )
    _spec = _importlib.spec_from_file_location(
        "uniteractive_update", _script_dir
    )
    if _spec is None or _spec.loader is None:
        raise ImportError("Submodule RNALoops was not correctly cloned.")
    _motifs = _importlib.module_from_spec(_spec)
    _spec.loader.exec_module(_motifs)
except ImportError as e:
    print(e)


_RNAMOTIFOLD_CONFIG_DIR: _Path = _Path.joinpath(
    RNAMOTIFOLD_ROOT_DIR, "configs"
)
RNAMOTIFOLD_DEFAULTS_CONFIG = _Path.joinpath(
    _RNAMOTIFOLD_CONFIG_DIR, "defaults.ini"
)
RNAMOTIFOLD_INSTALL_CONFIG: _Path = _Path.joinpath(
    _RNAMOTIFOLD_CONFIG_DIR, "install.ini"
)
RNAMOTIFOLD_BIN: _Path = _Path.joinpath(RNAMOTIFOLD_ROOT_DIR, "bin")
RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
RNALOOPS_PATH: _Path = _check_submodule("RNALoops")
_RNAMOTIFOLD_MOTIFS_PATH: _Path = _Path.joinpath(
    RNALOOPS_PATH,
    "Misc",
    "Applications",
    "RNAmotiFold",
    "motifs",
    "versions",
    "combined",
)
AVAILABLE_BINARIES: list[str] = [
    "RNAmotiFold",
    "RNAmoSh",
    "RNAmotiCes",
    "RNAmotiFoldMotmicro",
    "RNAmoShMotmicro",
    "RNAmotiCesMotmicro",
    "RNAmotiFold_pfc",
    "RNAmoSh_pfc",
    "RNAmotiCes_pfc",
    "RNAmotiFold_motmacro_pfc",
    "RNAmoSh_motmacro_pfc",
    "RNAmotiCes_motmacro_pfc",
    "RNAmotiFold_subopt",
    "RNAmoSh_subopt",
    "RNAmotiCes_subopt",
    "RNAmotiFold_motmacro_subopt",
    "RNAmoSh_motmacro_subopt",
    "RNAmotiCes_motmacro_subopt",
    "RNAmotiAlign",
]
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
