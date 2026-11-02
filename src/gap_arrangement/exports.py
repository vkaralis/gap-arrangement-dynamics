"""Export figures without embedded metadata and validate their actual content."""
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from PIL import Image


def save_figure(figure, destination):
    """Preserve plotted pixels/vectors while omitting authoring metadata."""
    destination = Path(destination)
    figure.savefig(destination.with_suffix('.pdf'), metadata={
        'Creator': None, 'Producer': None, 'CreationDate': None,
        'ModDate': None, 'Author': None, 'Title': None, 'Subject': None,
        'Keywords': None, 'Trapped': None,
    })
    png = destination.with_suffix('.png')
    figure.savefig(png, metadata={'Software': None})
    # New image from pixels: no text, EXIF, resolution, ICC or software fields.
    # The numerical figure is not resized or recolored.
    with Image.open(png) as source:
        pixels = source.convert('RGBA').tobytes()
        size = source.size
    Image.frombytes('RGBA', size, pixels).save(png)


def check_csv(path, fields, expected):
    """Check column names/order, row count and every numeric cell against JSON."""
    path = Path(path)
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == fields, f'{path.name}: incorrect CSV columns'
        rows = list(reader)
    assert len(rows) == len(expected), f'{path.name}: incorrect CSV row count'
    for index, (actual, wanted) in enumerate(zip(rows, expected)):
        assert set(actual) == set(fields), f'{path.name}: malformed row {index + 1}'
        for key in fields:
            try:
                value = float(actual[key])
            except (TypeError, ValueError) as error:
                raise AssertionError(f'{path.name}: nonnumeric {key} in row {index + 1}') from error
            assert np.isfinite(value), f'{path.name}: non-finite {key}'
            assert np.isclose(value, wanted[key], rtol=1e-12, atol=0), (
                f'{path.name}: changed {key} in row {index + 1}')


def check_figures(actual_dir, expected_dir):
    """Compare freshly rendered content in the current software environment."""
    actual_dir, expected_dir = Path(actual_dir), Path(expected_dir)
    for name in ['figure1_approximation', 'figure2_separation']:
        with Image.open(actual_dir / f'{name}.png') as actual, Image.open(expected_dir / f'{name}.png') as expected:
            assert not actual.info, f'{name}.png: embedded metadata is present'
            assert actual.mode == expected.mode and actual.size == expected.size, (
                f'{name}.png: changed image mode or dimensions')
            np.testing.assert_array_equal(np.asarray(actual), np.asarray(expected),
                                          err_msg=f'{name}.png: plotted content differs')
        # Metadata-free PDFs are deterministic for these fixed plots within one
        # rendering environment. This detects swapped, truncated or edited PDFs.
        assert (actual_dir / f'{name}.pdf').read_bytes() == (expected_dir / f'{name}.pdf').read_bytes(), (
            f'{name}.pdf: exported content differs; regenerate with this environment')


def verify_exports(output_dir, result, data_dir):
    from .simulation import make_outputs

    output_dir, data_dir = Path(output_dir), Path(data_dir)
    for name, rows in [('table1_particles.csv', result['particle_rows']),
                       ('table2_response.csv', result['observable_rows'])]:
        check_csv(output_dir / 'tables' / name, list(rows[0]), rows)
    check_csv(output_dir / 'tables/theoretical_coefficients.csv', ['beta', 'C2_over_2'], [
        {'beta': 0., 'C2_over_2': result['affine_C2_over_2']},
        {'beta': result['parameters']['beta'], 'C2_over_2': result['C2_over_2']},
    ])
    with TemporaryDirectory(prefix='gap-export-check-') as temporary:
        with np.load(data_dir / 'simulation_curves.npz', allow_pickle=False) as curves:
            make_outputs(result, curves, data_dir, Path(temporary))
        check_figures(output_dir / 'figures', Path(temporary) / 'figures')
    return {'csv_values': 'passed', 'png_pixels': 'passed', 'pdf_content': 'passed',
            'png_metadata': 'absent', 'comparison': 'rerendered from validated data in the current environment'}
