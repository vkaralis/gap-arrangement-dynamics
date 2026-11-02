from contextlib import redirect_stderr
import io
from pathlib import Path
import re
import json
import tomllib
import unittest
from unittest.mock import patch

from gap_arrangement import __version__
from gap_arrangement.cli import main
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]


class ReleaseMetadataTests(unittest.TestCase):
    def test_versions_and_python_support_agree(self):
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
        self.assertEqual(project["version"], __version__)
        self.assertEqual(project["requires-python"], ">=3.12,<3.13")
        self.assertEqual((ROOT / ".python-version").read_text(encoding="utf-8").strip(), "3.12")
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        self.assertEqual(re.search(r"^version: (.+)$", citation, re.M).group(1), __version__)
        self.assertIn(f"## [{__version__}]", (ROOT / "CHANGELOG.md").read_text(encoding="utf-8"))
        for workflow in (ROOT / ".github/workflows").glob("*.yml"):
            self.assertIn("python-version-file: '.python-version'", workflow.read_text(encoding="utf-8"))

    def test_nine_topics_agree_with_package_keywords(self):
        topics = (ROOT / "topics.txt").read_text(encoding="utf-8").splitlines()
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
        self.assertEqual(len(topics), 9)
        self.assertEqual(topics, project["keywords"])
        citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
        for topic in topics:
            self.assertIn(f"  - {topic}\n", citation)

    def test_source_entry_point_rejects_unsupported_python(self):
        error = io.StringIO()
        with patch("gap_arrangement.cli.sys.version_info", (3, 13, 0)), redirect_stderr(error):
            with self.assertRaises(SystemExit) as caught:
                main(["certificate"])
        self.assertEqual(caught.exception.code, 2)
        self.assertIn("supports Python 3.12 only", error.getvalue())

    def test_code_creator_and_study_contributor_are_distinct(self):
        project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
        self.assertEqual(project['authors'], [{'name': 'Vangelis D. Karalis'}])
        citation = (ROOT / 'CITATION.cff').read_text(encoding="utf-8")
        self.assertIn('family-names: "Karalis"', citation)
        self.assertIn('given-names: "Vangelis D."', citation)
        self.assertNotIn('Georgia', citation)
        readme = (ROOT / 'README.md').read_text(encoding="utf-8")
        self.assertIn('**Code creator:** Vangelis D. Karalis', readme)
        self.assertIn('**Study contributor:** Georgia Karali', readme)
        self.assertIn('Department of Mathematics, National and\nKapodistrian University of Athens', readme)

    def test_public_reference_metadata_removed(self):
        for path in (ROOT / 'reference/data').glob('*.json'):
            record = json.loads(path.read_text(encoding="utf-8"))
            self.assertNotIn('environment', record)
            self.assertNotIn('elapsed_seconds', record)
        for path in (ROOT / 'reference/figures').glob('*.png'):
            with Image.open(path) as image:
                self.assertFalse(image.info, path.name)
        for path in (ROOT / 'reference/figures').glob('*.pdf'):
            for field in [b'/Creator', b'/Producer', b'/CreationDate', b'/ModDate', b'/Author', b'/Metadata']:
                self.assertNotIn(field, path.read_bytes(), path.name)
