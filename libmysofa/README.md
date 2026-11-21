# vendored libmysofa

Trimmed copy of libmysofa (https://github.com/hoene/libmysofa), used to parse
AES69 SOFA files at runtime. Only the HDF5 container reader and the pieces it
needs are included; the resampling, interpolation, and k-d tree lookup code is
not.

- Upstream revision: 648eed03472e6720a1ea45d1a1f86c4efb569ff9 (v1.3.5)
- License: BSD-3-Clause, see `LICENSE`.

## Layout

- `hdf/` — the HDF5 reader (superblock, object headers, fractal heap, B-tree,
  global heap, deflate via zlib).
- `hrtf/reader.c` — builds the `MYSOFA_HRTF` from the HDF5 datasets.
- `hrtf/check.c` — SOFA convention validation (used for a warning only).
- `hrtf/spherical.c`, `hrtf/tools.c` — cartesian/spherical coordinate helpers.
- `hrtf/config.h`, `hrtf/mysofa_export.h` — minimal replacements for the files
  upstream generates with CMake; `MYSOFA_EXPORT` expands to nothing since the
  library is statically linked.

## Local changes vs upstream

- Include paths adjusted for the flattened layout.
- `hdf/dataobject.c`: `readOHDRHeaderMessageContinue` now decrements
  `recursive_counter` after recursing, so the counter tracks recursion depth
  instead of the total number of object-header continuations. Recent
  netCDF-4/HDF5 writers emit enough continuations that valid files otherwise
  hit the limit and fail with `MYSOFA_UNSUPPORTED_FORMAT`.

The code is built as the static library `alsoft.mysofa` with warnings off (it
is third-party C code written against a different warning set). It is linked
into the main library only when `ALSOFT_RUNTIME_SOFA` is enabled and zlib is
available.
