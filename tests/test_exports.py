import csv
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import numpy as np
from PIL import Image, PngImagePlugin

from gap_arrangement.exports import check_csv, check_figures, verify_exports
from gap_arrangement.simulation import make_outputs

ROOT = Path(__file__).resolve().parents[1]


class ExportContentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = tempfile.TemporaryDirectory()
        cls.expected = Path(cls.workspace.name) / 'expected'
        cls.data = ROOT / 'reference/data'
        cls.result = json.loads((cls.data / 'simulation_results.json').read_text())
        with np.load(cls.data / 'simulation_curves.npz', allow_pickle=False) as arrays:
            make_outputs(cls.result, arrays, cls.data, cls.expected)

    @classmethod
    def tearDownClass(cls):
        cls.workspace.cleanup()

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.output = Path(self.temporary.name) / 'outputs'
        shutil.copytree(self.expected, self.output)

    def tearDown(self):
        self.temporary.cleanup()

    def test_correct_exports_pass_full_content_check(self):
        self.assertEqual(verify_exports(self.output, self.result, self.data)['pdf_content'], 'passed')

    def test_csv_changed_value_rejected(self):
        path = self.output / 'tables/table1_particles.csv'
        with path.open(newline='') as handle:
            reader = csv.DictReader(handle)
            fields, rows = reader.fieldnames, list(reader)
        rows[0]['max_mesh_E_plus_D'] = '999'
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader(); writer.writerows(rows)
        with self.assertRaisesRegex(AssertionError, 'changed max_mesh_E_plus_D'):
            verify_exports(self.output, self.result, self.data)

    def test_csv_missing_row_rejected(self):
        path = self.output / 'tables/table2_response.csv'
        path.write_text('\n'.join(path.read_text().splitlines()[:-1]) + '\n')
        with self.assertRaisesRegex(AssertionError, 'row count'):
            verify_exports(self.output, self.result, self.data)

    def test_csv_changed_header_rejected(self):
        path = self.output / 'tables/theoretical_coefficients.csv'
        path.write_text(path.read_text().replace('C2_over_2', 'wrong_column'))
        with self.assertRaisesRegex(AssertionError, 'columns'):
            verify_exports(self.output, self.result, self.data)

    def test_csv_nan_rejected(self):
        path = self.output / 'tables/theoretical_coefficients.csv'
        path.write_text('beta,C2_over_2\n0,nan\n0.5,0.741984411175\n')
        with self.assertRaisesRegex(AssertionError, 'non-finite'):
            verify_exports(self.output, self.result, self.data)

    def test_swapped_png_rejected(self):
        figures = self.output / 'figures'
        shutil.copyfile(figures / 'figure2_separation.png', figures / 'figure1_approximation.png')
        with self.assertRaises(AssertionError):
            check_figures(figures, self.expected / 'figures')

    def test_modified_png_pixel_rejected(self):
        path = self.output / 'figures/figure1_approximation.png'
        with Image.open(path) as source:
            image = source.copy()
        image.putpixel((0, 0), (0, 0, 0, 255))
        image.save(path)
        with self.assertRaises(AssertionError):
            check_figures(self.output / 'figures', self.expected / 'figures')

    def test_swapped_pdf_rejected(self):
        figures = self.output / 'figures'
        shutil.copyfile(figures / 'figure2_separation.pdf', figures / 'figure1_approximation.pdf')
        with self.assertRaisesRegex(AssertionError, 'figure1_approximation.pdf'):
            check_figures(figures, self.expected / 'figures')

    def test_png_metadata_rejected(self):
        path = self.output / 'figures/figure1_approximation.png'
        with Image.open(path) as source:
            image = source.copy()
        info = PngImagePlugin.PngInfo(); info.add_text('Author', 'test metadata')
        image.save(path, pnginfo=info)
        with self.assertRaisesRegex(AssertionError, 'metadata'):
            check_figures(self.output / 'figures', self.expected / 'figures')
