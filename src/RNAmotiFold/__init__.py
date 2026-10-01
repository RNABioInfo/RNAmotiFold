from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location

_RNAMOTIFOLD_ROOT_DIR = Path(__file__).resolve().absolute().parent

try:
    script_dir = (
        _RNAMOTIFOLD_ROOT_DIR
        / "RNALoops"
        / "Misc"
        / "Applications"
        / "RNAmotiFold"
        / "motifs"
        / "get_RNA3D_motifs.py"
    )
    spec = spec_from_file_location("uniteractive_update", script_dir)
    if spec is None or spec.loader is None:
        raise ImportError(
            f"Submodule RNALoops was not correctly cloned."
        )
    motifs = module_from_spec(spec)
    spec.loader.exec_module(motifs)
except ImportError as e:
    raise e