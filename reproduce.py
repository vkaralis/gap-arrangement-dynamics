"""Source-checkout entry point; run `python reproduce.py --help` for commands."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from gap_arrangement.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
