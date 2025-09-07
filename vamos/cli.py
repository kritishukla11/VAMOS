import argparse
from .utils import setup_logging
from .mapping import run_mapping
from .clustering import run_clustering
from .ml import run_ml
from .stats import run_stats

def main():
    parser = argparse.ArgumentParser(
        prog="vamos",
        description="VAMOS: Variant Mapping and Oncogenic Signatures toolkit"
    )
    parser.add_argument("-v", "--verbose", action="count", default=1, help="Increase verbosity (-v, -vv)")
parser.add_argument("-c", "--config", help="Path to config.yaml", default=None)


    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("map", help="Run AlphaFold mapping step")
    subparsers.add_parser("cluster", help="Run density-based clustering")
    subparsers.add_parser("ml", help="Run machine learning pipeline")
    subparsers.add_parser("stats", help="Run statistical analyses")

    args, unknown = parser.parse_known_args()
    from .config import Config
setup_logging(args.verbose)
cfg = Config(args.config) if args.config else Config()


    if args.command == "map":
        run_mapping()
    elif args.command == "cluster":
        run_clustering()
    elif args.command == "ml":
        run_ml()
    elif args.command == "stats":
        run_stats()

if __name__ == "__main__":
    main()
