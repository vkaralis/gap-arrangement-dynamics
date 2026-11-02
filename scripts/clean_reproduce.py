"""Build and reproduce in a new isolated Python 3.12 virtual environment."""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import venv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, default=Path("clean-run"),
                        help="A new directory; existing directories are refused")
    args = parser.parse_args()
    if sys.version_info[:2] != (3, 12):
        parser.error("Use Python 3.12 for this release.")
    source = Path(__file__).resolve().parents[1]
    work = args.work_dir.resolve()
    if work.exists():
        parser.error("Choose a new --work-dir; an existing environment is not a clean test.")
    start = time.monotonic()
    started = datetime.now(timezone.utc).isoformat()
    work.mkdir(parents=True)
    checkout = work / "source"

    def ignored(directory, names):
        excluded = {".git", ".venv", "venv", "__pycache__", "build", "dist", "outputs",
                    "clean-run", ".pytest_cache", ".DS_Store", "generated", "input-check"}
        return [name for name in names if name in excluded or name.endswith((".pyc", ".egg-info"))
                or name.startswith("clean-run-")
                or (Path(directory) / name).resolve() == work]

    shutil.copytree(source, checkout, ignore=ignored)
    environment = os.environ.copy()
    for key in ["PYTHONPATH", "PYTHONHOME", "PYTHONOPTIMIZE", "VIRTUAL_ENV"]:
        environment.pop(key, None)
    environment["PYTHONNOUSERSITE"] = "1"
    environment["PIP_DISABLE_PIP_VERSION_CHECK"] = "1"
    virtual = work / "venv"
    print("Creating a new venv with system-site-packages disabled", flush=True)
    venv.EnvBuilder(with_pip=True, system_site_packages=False).create(virtual)
    python = virtual / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    outside = work / "outside"
    outside.mkdir()

    def run(arguments, log_name):
        print(f"Running {log_name}", flush=True)
        with (work / log_name).open("w", encoding="utf-8") as log:
            process = subprocess.run([str(python), "-I"] + arguments, cwd=outside,
                                     env=environment, stdout=log, stderr=subprocess.STDOUT)
        if process.returncode:
            raise RuntimeError(f"Command failed; see {work / log_name}")

    try:
        run(["-m", "pip", "--isolated", "install", "--no-cache-dir", "--index-url",
             "https://pypi.org/simple", "-r", str(checkout / "requirements-lock.txt"),
             "-r", str(checkout / "requirements-build.txt"), "--report",
             str(work / "dependency_install.json")], "install.log")
        run(["-m", "pip", "--isolated", "install", "--no-index", "--no-cache-dir",
             "--no-build-isolation", "--no-deps", str(checkout)], "build.log")
        run(["-m", "pip", "check"], "pip-check.log")
        proof = r'''
from importlib.metadata import distributions
from pathlib import Path
import json, site, sys
import numpy, scipy, matplotlib, gap_arrangement
prefix = Path(sys.prefix).resolve()
assert sys.prefix != sys.base_prefix
assert site.ENABLE_USER_SITE is False
assert 'include-system-site-packages = false' in (prefix/'pyvenv.cfg').read_text()
modules = [numpy, scipy, matplotlib, gap_arrangement]
assert all(prefix in Path(module.__file__).resolve().parents for module in modules)
assert all(prefix in Path(path).resolve().parents or Path(path).resolve() == prefix
           for path in sys.path if 'site-packages' in path)
installed = {d.metadata['Name']: d.version for d in distributions()}
assert all(prefix in Path(d.locate_file('')).resolve().parents for d in distributions())
print(json.dumps({'system_site_packages': False, 'user_site_enabled': site.ENABLE_USER_SITE,
 'all_distributions_inside_venv': True, 'all_numerical_modules_inside_venv': True,
 'python': sys.version, 'software_version': gap_arrangement.__version__,
 'packages': installed}, indent=2))
'''
        run(["-c", proof], "isolation.json")
        run(["-m", "unittest", "discover", "-s", str(checkout / "tests"), "-v"], "tests.log")
        run(["-m", "gap_arrangement", "reproduce", "--reference-dir", str(checkout / "reference"),
             "--inputs-dir", str(checkout / "inputs"), "--output-dir", str(work / "results")],
            "reproduction.log")
        manifest = json.loads((work / "results/run_manifest.json").read_text())
        assert manifest["status"] == "passed"
        for name, expected in manifest["source_sha256"].items():
            assert sha256((checkout / "src/gap_arrangement" / name).read_bytes()).hexdigest() == expected
        tests = (work / "tests.log").read_text()
        count = int(re.search(r"Ran (\d+) tests?", tests).group(1))
        result = {"status": "passed", "started_utc": started,
                  "completed_utc": datetime.now(timezone.utc).isoformat(),
                  "elapsed_seconds": time.monotonic() - start,
                  "isolation": json.loads((work / "isolation.json").read_text()),
                  "tests_passed": count,
                  "input_verification": json.loads((work / "results/data/input_verification.json").read_text()),
                  "numerical_verification": json.loads((work / "results/data/verification.json").read_text()),
                  "source_sha256": manifest["source_sha256"],
                  "dependency_source": "Fresh downloads from https://pypi.org/simple; no pip cache",
                  "execution": "Installed wheel, isolated Python (-I), working directory outside source"}
    except (Exception, KeyboardInterrupt) as error:
        (work / "clean_reproduction.json").write_text(json.dumps({"status": "failed",
            "started_utc": started, "error": str(error)}, indent=2) + "\n", encoding="utf-8")
        raise
    (work / "clean_reproduction.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Clean reproduction passed: {work / 'clean_reproduction.json'}", flush=True)


if __name__ == "__main__":
    main()
