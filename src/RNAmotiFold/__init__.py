from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location

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
    script_dir: Path = (
        RNAMOTIFOLD_ROOT_DIR
        / "RNALoops"
        / "Misc"
        / "Applications"
        / "RNAmotiFold"
        / "motifs"
        / "get_RNA3D_motifs.py"
    )
    spec = spec_from_file_location("uniteractive_update", script_dir)
    if spec is None or spec.loader is None:
        raise ImportError(f"Submodule RNALoops was not correctly cloned.")
    motifs = module_from_spec(spec)
    spec.loader.exec_module(motifs)
except ImportError as e:
    raise e


RNAMOTIFOLD_CONFIG_DIR: Path = Path.joinpath(RNAMOTIFOLD_ROOT_DIR, "configs")
RNAMOTIFOLD_DEFAULTS_CONFIG = Path.joinpath(RNAMOTIFOLD_CONFIG_DIR, "defaults.ini")
RNAMOTIFOLD_PATHS_CONFIG: Path = Path.joinpath(RNAMOTIFOLD_CONFIG_DIR, "paths.ini")
RNAMOTIFOLD_BIN: Path = Path.joinpath(RNAMOTIFOLD_ROOT_DIR, "bin")
RNAMOTIFOLD_BIN.mkdir(exist_ok=True, parents=True)
RNALOOPS_PATH: Path = _check_submodule("RNALoops")
RNAMOTIFOLD_MOTIFS_PATH: Path = Path.joinpath(
    RNALOOPS_PATH, "Misc", "Applications", "RNAmotiFold", "motifs", "versions", "combined"
)
AVAILABLE_ALGORITHMS: list[str] = [
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
AVAILABLE_VERSIONS: list[str] = [x.name for x in RNAMOTIFOLD_MOTIFS_PATH.iterdir() if x.is_dir()]
