"""Command-line interface for the fixed numerical experiments."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

from . import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", nargs="?", default="reproduce", choices=[
        "reproduce", "validate", "simulate", "affine", "certificate", "figures", "verify",
        "synthetic", "check-inputs"])
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"),
                        help="Generated files (default: ./outputs)")
    parser.add_argument("--reference-dir", type=Path, default=Path("reference"),
                        help="Frozen reference directory, needed by reproduce/verify")
    parser.add_argument("--inputs-dir", type=Path, default=Path("inputs"),
                        help="Bundled synthetic data, checked by reproduce/check-inputs")
    parser.add_argument("--input-dir", type=Path,
                        help="For figures only: input directory containing data/ (default: output directory)")
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args(argv)
    if sys.version_info[:2] != (3, 12):
        parser.error("This release supports Python 3.12 only (>=3.12,<3.13).")
    if sys.flags.optimize:
        parser.error("Run without -O or PYTHONOPTIMIZE: numerical assertions must remain enabled.")
    if args.input_dir is not None and args.command != "figures":
        parser.error("--input-dir is used only by the figures command.")
    output, reference = args.output_dir.resolve(), args.reference_dir.resolve()
    source_reference = Path(__file__).resolve().parents[2] / "reference"
    source_inputs = Path(__file__).resolve().parents[2] / "inputs"
    for protected in [reference, source_reference.resolve(), args.inputs_dir.resolve(), source_inputs.resolve()]:
        if output == protected or protected in output.parents:
            parser.error("The output directory cannot overwrite frozen reference or input files.")
    if args.command in ["reproduce", "verify"] and not (reference / "data/simulation_results.json").is_file():
        parser.error("Reference data not found. Pass --reference-dir /path/to/repository/reference.")
    started = datetime.now(timezone.utc).isoformat()
    output.mkdir(parents=True, exist_ok=True)
    # A previous successful manifest must not survive an interrupted new command.
    manifest = output / "run_manifest.json"
    manifest.write_text(json.dumps({"command": args.command, "status": "running",
                                   "started_utc": started}, indent=2) + "\n", encoding="utf-8")
    try:
        if args.command in ["reproduce", "check-inputs"]:
            from .synthetic import validate_inputs
            data = output / "data"
            data.mkdir(parents=True, exist_ok=True)
            result = validate_inputs(args.inputs_dir)
            (data / "input_verification.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            print(json.dumps(result, indent=2), flush=True)
        if args.command == "synthetic":
            from .synthetic import write_inputs
            print(json.dumps(write_inputs(output / "inputs"), indent=2), flush=True)
        if args.command in ["reproduce", "validate"]:
            from .validation import main as validate
            validate(output)
        if args.command in ["reproduce", "simulate"]:
            from .simulation import main as simulate
            simulate(output)
        if args.command in ["reproduce", "affine"]:
            from .affine import main as affine
            affine(output)
        if args.command in ["reproduce", "certificate"]:
            from .certificate import certificate
            data = output / "data"
            data.mkdir(parents=True, exist_ok=True)
            result = certificate()
            (data / "certificate.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            print(json.dumps(result, indent=2), flush=True)
        if args.command == "figures":
            from .simulation import make_outputs
            import numpy as np
            data = (args.input_dir.resolve() if args.input_dir else output) / "data"
            result = json.loads((data / "simulation_results.json").read_text())
            with np.load(data / "simulation_curves.npz", allow_pickle=False) as curves:
                make_outputs(result, curves, data, output)
        if args.command in ["reproduce", "verify"]:
            from .reproducibility import verify
            print(json.dumps(verify(output, reference), indent=2), flush=True)
        # The rational certificate command has no third-party dependencies.
        if args.command == "certificate":
            manifest.write_text(json.dumps({"command": "certificate", "status": "passed",
                "started_utc": started, "completed_utc": datetime.now(timezone.utc).isoformat(),
                "python": sys.version}, indent=2) + "\n", encoding="utf-8")
        else:
            from .reproducibility import write_manifest
            write_manifest(output, args.command, started, reference, args.inputs_dir)
        print(f"Completed {args.command}: {output}", flush=True)
        return 0
    except (Exception, KeyboardInterrupt) as error:
        manifest.write_text(json.dumps({"command": args.command, "status": "failed",
            "started_utc": started, "error": str(error)}, indent=2) + "\n", encoding="utf-8")
        raise
