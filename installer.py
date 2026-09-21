import src.RNAmotiFold.bgap_rna.alg_setup
from pathlib import Path
try:
    import submodules.RNALoops.Misc.Applications.RNAmotiFold.motifs.get_RNA3D_motifs as motifs
except ImportError as e:
    print(
        f"Submodule was not correctly cloned. If you didn't clone this repo with --recurse-submodules run git submodule update --init --recursive from {Path(__file__).parent.absolute()}"
    )
    raise e

def main():
    """main setup function that checks for the gap compiler, installs it if necessary, fetches newest motif sequences and (re)compiles all preset algorithms (RNAmotiFold, RNAmoSh, RNAmotiCes)"""
    args = src.RNAmotiFold.bgap_rna.alg_setup.get_cmd_args()
    done: bool = False
    motifs.uninteractive_update(args.version)  # type: ignore

    done = src.RNAmotiFold.bgap_rna.alg_setup.setup_algorithms(
        args.gapc_path, args.perl_path, int(args.workers)
    )
    if done:
        print("Algorithms are all set up, you can now use RNAmotiFold")
    else:
        print("Something went wrong compiling the RNAmotiFold algorithms, please check outputs")

if __name__ == "__main__":
    main()
