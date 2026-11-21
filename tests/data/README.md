# SOFA test fixtures

Small AES69 SOFA files used by `tests/sofa_hrtf_test` (and by the loader
tests it drives). They are committed so the test doesn't need Python, h5py, or
netCDF4 to run.

`gen_sofa_fixtures.py` regenerates them. It uses netCDF4 to write netCDF-4
(HDF5) files, which is what real SOFA files use.

```sh
# Requires Python with netCDF4 and numpy. The HDF5 version linked into
# netCDF4 must be one the bundled libmysofa reader understands (1.8-1.14);
# HDF5 2.x writes structures it doesn't parse. netCDF4 1.6.x bundles HDF5
# 1.12, which works.
python3 gen_sofa_fixtures.py .
```

Contents:

- `valid.sofa` — stereo, 44.1kHz, 8 azimuths x 9 elevations, I,R delays.
- `valid_mono.sofa` — mono (R=1), 48kHz, no delay data (onset is used).
- `valid_mr_delay.sofa` — stereo, M,R delays (0.0001s left, 0.0002s right).
- `long_ir.sofa` — 256-sample responses (> `HrirLength`), exercises windowing.
- `bad_receivers.sofa` — R=3, rejected as a dimension error.
- `bad_rate.sofa` — 1kHz sample rate, rejected as a dimension error.
- `bad_grid.sofa` — too few elevations to form a uniform grid.
- `missing_ir.sofa` — no `Data.IR`, missing-data error.
- `missing_sr.sofa` — no `Data.SamplingRate`, missing-data error.
- `silent_ir.sofa` — all-zero responses, missing-data error.
- `truncated.sofa` — HDF5 signature followed by garbage.
- `not_sofa.mhr` — an `.mhr`-like marker with a truncated body; goes through
  the native loader and must not raise a SOFA error.
