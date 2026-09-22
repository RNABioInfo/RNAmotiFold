import src.RNAmotiFold.bgap_rna.alg_setup
from pathlib import Path

try:
    import submodules.RNALoops.Misc.Applications.RNAmotiFold.motifs.get_RNA3D_motifs as motifs
except ImportError as e:
    print(
        f"Submodules were not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {Path(__file__).parent.absolute()}"
    )
    raise e

if __name__ == "__main__":
    src.RNAmotiFold.bgap_rna.alg_setup.main()
