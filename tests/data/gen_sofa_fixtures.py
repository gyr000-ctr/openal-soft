#!/usr/bin/env python3
"""Generates the AES69 SOFA test fixtures used by sofa_hrtf_test.

The fixtures are netCDF-4 (HDF5) files, written through the netCDF4 library so
the netCDF dimension scales and DIMENSION_LIST attributes match what real SOFA
files contain. They are small (a handful of elevations and azimuths) and are
committed so the test doesn't need Python to run.

Run with: python3 gen_sofa_fixtures.py <output-dir>

Requires netCDF4 and numpy. The HDF5 version linked into netCDF4 must be old
enough for the bundled libmysofa reader (it handles HDF5 1.8-1.14; HDF5 2.x
files use structures it doesn't understand).
"""

import os
import sys

import netCDF4
import numpy as np


def grid_positions(az_count, ev_count, radius=1.4):
    """Returns (M, 3) spherical source positions (azimuth 0 = front, increasing
    toward the left; elevation -90..90; distance in metres), matching the
    OpenAL Soft polar layout. Poles are single points.
    """
    rows = []
    for ei in range(ev_count):
        elev = -90.0 + 180.0 * ei / (ev_count - 1)
        if abs(elev) >= 89.999:
            rows.append((0.0, elev, radius))
        else:
            for ai in range(az_count):
                rows.append((360.0 * ai / az_count, elev, radius))
    return np.array(rows, dtype=np.float64)


def write_sofa(path, positions, n=32, rate=44100.0, r=2, ir=None,
               delay=None, delay_type="I_R", extra_attrs=None,
               skip_ir=False, skip_rate=False):
    m = positions.shape[0]
    nc = netCDF4.Dataset(path, "w", format="NETCDF4")

    attrs = {
        "Conventions": "SOFA",
        "Version": "1.0",
        "SOFAConventions": "SimpleFreeFieldHRIR",
        "SOFAConventionsVersion": "1.0",
        "DataType": "FIR",
        "RoomType": "free field",
        "Title": "Runtime SOFA loader test fixture",
    }
    if extra_attrs:
        attrs.update(extra_attrs)
    for key, value in attrs.items():
        if value is not None:
            nc.setncattr(key, value)

    nc.createDimension("I", 1)
    nc.createDimension("C", 3)
    nc.createDimension("E", 1)
    nc.createDimension("M", m)
    nc.createDimension("R", r)
    nc.createDimension("N", n)

    def var(name, dims, data, vattrs=None):
        v = nc.createVariable(name, "f8", dims)
        v[:] = data
        if vattrs:
            for key, value in vattrs.items():
                v.setncattr(key, value)
        return v

    var("ListenerPosition", ("I", "C"), np.zeros((1, 3)), {"Type": "cartesian"})
    var("ListenerUp", ("I", "C"), np.array([[0.0, 0.0, 1.0]]), {"Type": "cartesian"})
    var("ListenerView", ("I", "C"), np.array([[1.0, 0.0, 0.0]]), {"Type": "cartesian"})
    var("ReceiverPosition", ("R", "C", "I"), np.zeros((r, 3, 1)),
        {"Type": "cartesian", "Units": "metre"})
    var("EmitterPosition", ("E", "C", "I"), np.zeros((1, 3, 1)),
        {"Type": "cartesian", "Units": "metre"})
    var("SourcePosition", ("M", "C"), positions,
        {"Type": "spherical", "Units": "degree, degree, metre"})

    if not skip_rate:
        var("Data.SamplingRate", ("I",), np.array([rate]), {"Units": "hertz"})

    if not skip_ir:
        if ir is None:
            ir = np.zeros((m, r, n), dtype=np.float64)
            for si in range(m):
                for ch in range(r):
                    # A short pre-ring and a clear peak at sample 8.
                    ir[si, ch, :8] = 0.01 * (si + 1) * (ch + 1)
                    ir[si, ch, 8] = 0.25 + 0.001 * si + 0.01 * ch
        var("Data.IR", ("M", "R", "N"), ir)

    if delay_type is not None and delay is not None:
        if delay_type == "I_R":
            var("Data.Delay", ("I", "R"), np.asarray(delay, dtype=np.float64))
        elif delay_type == "M_R":
            var("Data.Delay", ("M", "R"), np.asarray(delay, dtype=np.float64))

    nc.close()


def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[1]} <output-dir>")
    outdir = sys.argv[1]
    os.makedirs(outdir, exist_ok=True)
    j = os.path.join

    full = grid_positions(8, 9)  # 1 + 7*8 + 1 = 58 measurements.

    # Valid stereo, 44.1kHz, I,R delays of zero.
    write_sofa(j(outdir, "valid.sofa"), full, n=32, rate=44100.0, r=2,
               delay=np.zeros((1, 2)), delay_type="I_R")

    # Valid mono (left ear only), 48kHz, no delay data at all.
    write_sofa(j(outdir, "valid_mono.sofa"), full, n=32, rate=48000.0, r=1,
               delay=None, delay_type=None)

    # Valid stereo with per-measurement delays (M,R). 0.0001s and 0.0002s at
    # 44.1kHz become ~4.4 and ~8.8 samples, within the delay range.
    m = full.shape[0]
    mr_delays = np.zeros((m, 2), dtype=np.float64)
    mr_delays[:, 0] = 0.0001
    mr_delays[:, 1] = 0.0002
    write_sofa(j(outdir, "valid_mr_delay.sofa"), full, n=32, rate=44100.0, r=2,
               delay=mr_delays, delay_type="M_R")

    # Over-long responses (N=256 > 128), to exercise the windowing path.
    long_ir = np.zeros((m, 2, 256), dtype=np.float64)
    for si in range(m):
        for ch in range(2):
            long_ir[si, ch, 16:24] = 0.01
            long_ir[si, ch, 64] = 0.2 + 0.01 * ch
            long_ir[si, ch, 128] = 0.05  # Would be cut off by windowing.
    write_sofa(j(outdir, "long_ir.sofa"), full, n=256, rate=44100.0, r=2,
               ir=long_ir, delay=np.zeros((1, 2)), delay_type="I_R")

    # Unsupported receiver count.
    write_sofa(j(outdir, "bad_receivers.sofa"), full, n=32, rate=44100.0, r=3,
               delay=np.zeros((1, 3)), delay_type="I_R")

    # Sample rate outside the supported range.
    write_sofa(j(outdir, "bad_rate.sofa"), full, n=32, rate=1000.0, r=2,
               delay=np.zeros((1, 2)), delay_type="I_R")

    # Too few elevation rows to build a uniform grid.
    sparse = grid_positions(8, 2)
    write_sofa(j(outdir, "bad_grid.sofa"), sparse, n=32, rate=44100.0, r=2,
               delay=np.zeros((1, 2)), delay_type="I_R")

    # Missing the FIR data.
    write_sofa(j(outdir, "missing_ir.sofa"), full, n=32, rate=44100.0, r=2,
               skip_ir=True, delay=np.zeros((1, 2)), delay_type="I_R")

    # Missing the sample rate.
    write_sofa(j(outdir, "missing_sr.sofa"), full, n=32, rate=44100.0, r=2,
               skip_rate=True, delay=np.zeros((1, 2)), delay_type="I_R")

    # Present but all-silent FIR data.
    silent = np.zeros((m, 2, 32), dtype=np.float64)
    write_sofa(j(outdir, "silent_ir.sofa"), full, n=32, rate=44100.0, r=2,
               ir=silent, delay=np.zeros((1, 2)), delay_type="I_R")

    # An HDF5 signature followed by garbage.
    with open(j(outdir, "truncated.sofa"), "wb") as f:
        f.write(b"\x89HDF\r\n\x1a\n")
        f.write(b"\x00" * 16)

    # Looks like a native .mhr marker, with a truncated body.
    with open(j(outdir, "not_sofa.mhr"), "wb") as f:
        f.write(b"MinPHR03")
        f.write(b"\x00" * 16)

    print("wrote fixtures to", outdir)


if __name__ == "__main__":
    main()
