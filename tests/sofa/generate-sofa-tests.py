#!/usr/bin/env python3
"""Generates the synthetic SOFA test assets in this directory.

The files are small, fully synthetic SimpleFreeFieldHRIR data sets used by
tests/sofa.t.cpp to exercise OpenAL Soft's runtime SOFA HRTF loading:

  stereo-44100.sofa / stereo-48000.sofa / stereo-96000.sofa
      Stereo (R=2) sets at various sample rates. The 96kHz set uses longer
      IRs (N=160) to exercise onset windowing, the others use N=64.
  stereo-48000-sofa2.sofa
      Same as stereo-48000.sofa but tagged with SOFA conventions 2.0.
  mono-48000.sofa
      Mono (R=1) set, mirrored to the right ear at load time.
  mr-delay.sofa
      Stereo set with per-measurement (M,R) delays and cartesian source
      positions.
  corrupt-longir.sofa
      Valid structure but N exceeds the runtime IR length limit.
  corrupt-nonfinite.sofa / corrupt-silent.sofa / corrupt-huge.sofa
      stereo-48000.sofa with NaN samples, all-zero IRs, and absurd
      amplitudes respectively.
  corrupt-truncated.sofa
      The first 30% of stereo-48000.sofa's bytes.
  corrupt-not-sofa.sofa
      Not an HDF5 file at all.

Requires the netCDF4 and numpy Python packages (pip install netCDF4 numpy).
SOFA files are netCDF4 files; writing them with the netCDF library keeps the
HDF5 object structure compatible with MySOFA's minimal reader. Run with no
arguments; files are written next to this script.
"""

import os

import numpy as np
from netCDF4 import Dataset

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Elevation rows with measurements: -40 to +90 in 10 degree steps. The loader
# has to synthesize the missing -90..-50 rows.
ELEVATIONS = list(range(-40, 91, 10))
AZIMUTHS_PER_ROW = 24  # Uniform 15 degree steps on non-pole rows.
RADIUS = 1.5
MAX_ITD_SECONDS = 0.0005


def make_grid():
    """Returns (M,3) spherical source positions (az, el, radius) in degrees."""
    positions = []
    for el in ELEVATIONS:
        if abs(el) > 89.999:
            positions.append((0.0, float(el), RADIUS))
            continue
        for ai in range(AZIMUTHS_PER_ROW):
            positions.append((ai * 360.0 / AZIMUTHS_PER_ROW, float(el), RADIUS))
    return positions


def make_compact_grid():
    """A minimal grid (5 elevations x 8 azimuths) for oversized-IR tests."""
    positions = []
    for el in range(-90, 91, 45):
        if abs(el) > 89.999:
            positions.append((0.0, float(el), RADIUS))
            continue
        for ai in range(8):
            positions.append((ai * 45.0, float(el), RADIUS))
    return positions


def sph2cart(az_deg, el_deg, r):
    # Matches mysofa_s2c: x = cos(el)cos(az)r, y = cos(el)sin(az)r, z = sin(el)r.
    az = np.radians(az_deg)
    el = np.radians(el_deg)
    x = np.cos(el) * r
    return (x * np.cos(az), x * np.sin(az), np.sin(el) * r)


def make_hrir(seed, nsamples, onset, itd_samples):
    """Builds a stereo (2,N) HRIR: a decaying noise burst per ear, with the
    given base onset and interaural time difference (in samples)."""
    rng = np.random.default_rng(seed)
    hrir = np.zeros((2, nsamples), dtype=np.float64)
    for chan in range(2):
        # Left ear gets +itd/2, right ear -itd/2.
        start = int(round(onset + (itd_samples if chan == 0 else -itd_samples) / 2.0))
        start = max(start, 0)
        length = nsamples - start
        if length <= 0:
            continue
        noise = rng.standard_normal(length)
        envelope = np.exp(-np.arange(length) / (nsamples / 6.0))
        hrir[chan, start:] = noise * envelope
    # Normalize so the global peak is 0.8.
    peak = np.abs(hrir).max()
    if peak > 0:
        hrir *= 0.8 / peak
    return hrir


def write_sofa(path, rate, nsamples, channels, delay_dims="I,R", cartesian=False,
               sofa_version="1.0", ir_mutator=None, grid=None):
    positions = grid if grid is not None else make_grid()
    mcount = len(positions)

    with Dataset(path, "w", format="NETCDF4") as f:
        f.Conventions = "SOFA"
        f.Version = sofa_version
        f.SOFAConventions = "SimpleFreeFieldHRIR"
        f.SOFAConventionsVersion = sofa_version
        f.DataType = "FIR"
        f.RoomType = "free field"
        f.DatabaseName = "OpenAL Soft synthetic tests"
        f.Title = os.path.basename(path)
        f.APIName = "generate-sofa-tests.py"
        f.APIVersion = "1.0"

        f.createDimension("I", 1)
        f.createDimension("C", 3)
        f.createDimension("R", channels)
        f.createDimension("E", 1)
        f.createDimension("N", nsamples)
        f.createDimension("M", mcount)

        if cartesian:
            srcpos = np.array([sph2cart(*p) for p in positions], dtype=np.float64)
            src_type, src_units = "cartesian", "metre"
        else:
            srcpos = np.array(positions, dtype=np.float64)
            src_type, src_units = "spherical", "degree, degree, metre"

        var = f.createVariable("SourcePosition", "f8", ("M", "C"))
        var.Type = src_type
        var.Units = src_units
        var[:] = srcpos

        var = f.createVariable("ReceiverPosition", "f8", ("R", "C"))
        var.Type = "cartesian"
        var.Units = "metre"
        var[:] = np.zeros((channels, 3))

        var = f.createVariable("EmitterPosition", "f8", ("E", "C", "I"))
        var.Type = "cartesian"
        var.Units = "metre"
        var[:] = np.zeros((1, 3, 1))

        var = f.createVariable("ListenerPosition", "f8", ("I", "C"))
        var.Type = "cartesian"
        var.Units = "metre"
        var[:] = np.zeros((1, 3))

        var = f.createVariable("ListenerUp", "f8", ("I", "C"))
        var.Type = "cartesian"
        var.Units = "metre"
        var[:] = [[0.0, 0.0, 1.0]]

        var = f.createVariable("ListenerView", "f8", ("I", "C"))
        var.Type = "cartesian"
        var.Units = "metre"
        var[:] = [[1.0, 0.0, 0.0]]

        # Build the HRIRs. The interaural time difference varies with azimuth
        # like a simple spherical head.
        base_onset = max(4, nsamples // 16)
        # Keep per-measurement delays non-negative even at the largest ITD.
        delay_base = base_onset + MAX_ITD_SECONDS * rate / 2.0
        irs = np.zeros((mcount, channels, nsamples), dtype=np.float64)
        delays = np.zeros((mcount, channels), dtype=np.float64)
        for mi, (az, el, _radius) in enumerate(positions):
            itd = MAX_ITD_SECONDS * rate * np.sin(np.radians(az))
            seed = (mi * 7919 + int(rate) * 31 + nsamples) & 0xFFFFFFFF
            stereo = make_hrir(seed, nsamples, base_onset, itd)
            if delay_dims == "M,R":
                # Encode the ITD in Data.Delay instead of the IR position.
                stereo = np.roll(stereo, -base_onset, axis=1)
                delays[mi, 0] = (delay_base + itd / 2.0) / rate
                delays[mi, 1] = (delay_base - itd / 2.0) / rate
            for chan in range(channels):
                irs[mi, chan, :] = stereo[chan if channels > 1 else 0, :]

        if ir_mutator is not None:
            irs = ir_mutator(irs)

        f.createVariable("Data.IR", "f8", ("M", "R", "N"))[:] = irs

        var = f.createVariable("Data.SamplingRate", "f8", ("I",))
        var.Units = "hertz"
        var[:] = [float(rate)]

        if delay_dims == "M,R":
            f.createVariable("Data.Delay", "f8", ("M", "R"))[:] = delays
        else:
            f.createVariable("Data.Delay", "f8", ("I", "R"))[:] = np.zeros((1, channels))


def main():
    out = lambda name: os.path.join(OUT_DIR, name)

    write_sofa(out("stereo-44100.sofa"), 44100, 64, 2)
    write_sofa(out("stereo-48000.sofa"), 48000, 64, 2)
    write_sofa(out("stereo-96000.sofa"), 96000, 160, 2)
    write_sofa(out("stereo-48000-sofa2.sofa"), 48000, 64, 2, sofa_version="2.0")
    write_sofa(out("mono-48000.sofa"), 48000, 64, 1)
    write_sofa(out("mr-delay.sofa"), 48000, 64, 2, delay_dims="M,R", cartesian=True)

    # Structurally valid, but the IRs are beyond the runtime length limit.
    write_sofa(out("corrupt-longir.sofa"), 48000, 5000, 1, grid=make_compact_grid())

    def with_nan(irs):
        irs = irs.copy()
        irs[0, 0, 5] = np.nan
        return irs

    def with_zeros(irs):
        return np.zeros_like(irs)

    def with_huge(irs):
        return irs * 1.0e12

    write_sofa(out("corrupt-nonfinite.sofa"), 48000, 64, 2, ir_mutator=with_nan)
    write_sofa(out("corrupt-silent.sofa"), 48000, 64, 2, ir_mutator=with_zeros)
    write_sofa(out("corrupt-huge.sofa"), 48000, 64, 2, ir_mutator=with_huge)

    # Truncated file: keep the first 30% of a valid file's bytes.
    with open(out("stereo-48000.sofa"), "rb") as src:
        data = src.read()
    with open(out("corrupt-truncated.sofa"), "wb") as dst:
        dst.write(data[: len(data) * 3 // 10])

    with open(out("corrupt-not-sofa.sofa"), "w", encoding="ascii") as dst:
        dst.write("this is not a SOFA file\n")

    print("Generated SOFA test assets in", OUT_DIR)


if __name__ == "__main__":
    main()
