"""Exercise the user entry points from outside the checkout."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CommandLineTests(unittest.TestCase):
    def run_cli(self, args, cwd, optimize=False):
        env = os.environ.copy()
        env["PYTHONPATH"] = str(ROOT / "src")
        command = [sys.executable] + (["-O"] if optimize else [])
        return subprocess.run(command + ["-m", "gap_arrangement"] + args,
                              cwd=cwd, env=env, text=True, capture_output=True)

    def test_certificate_from_other_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            result = self.run_cli(["certificate"], temp)
            self.assertEqual(result.returncode, 0, result.stderr)
            data = json.loads((Path(temp) / "outputs/data/certificate.json").read_text())
            self.assertEqual(data["cap"], 10**9)

    def test_missing_audit_has_actionable_error_and_failed_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            result = self.run_cli(["simulate"], temp)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('run "validate" first', result.stderr)
            manifest = json.loads((Path(temp) / "outputs/run_manifest.json").read_text())
            self.assertEqual(manifest["status"], "failed")

    def test_reference_directory_cannot_be_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            result = self.run_cli(["certificate", "--output-dir", str(ROOT / "reference")], temp)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("cannot overwrite frozen reference", result.stderr)

    def test_optimized_python_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            result = self.run_cli(["certificate"], temp, optimize=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("assertions must remain enabled", result.stderr)


if __name__ == "__main__":
    unittest.main()
