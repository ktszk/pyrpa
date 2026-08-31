from ctypes import byref, c_int64, c_double
from pathlib import Path
import numpy as np

# Resolve lib path relative to this file so import works from any current working directory.
_LIB_DIR = Path(__file__).resolve().parents[1]


# -march target -> (rank, CPU flag it requires).  Higher rank = newer ISA, so the
# best library a CPU can run is the highest-ranked entry whose flag it has.  Only
# the targets `make auto` / `make dispatch` can emit need to be listed; an unknown
# suffix is treated as rank 0 with no requirement, i.e. usable but least preferred.
_ARCH_REQUIREMENTS = {
    "znver5":         (5, "avx512f"),
    "znver4":         (4, "avx512f"),
    "x86-64-v4":      (4, "avx512f"),
    "sapphirerapids": (4, "avx512f"),
    "icelake-server": (4, "avx512f"),
    "skylake-avx512": (4, "avx512f"),
    "znver3":         (2, "avx2"),
    "znver2":         (2, "avx2"),
    "znver1":         (2, "avx2"),
    "x86-64-v3":      (2, "avx2"),
    "haswell":        (2, "avx2"),
    "x86-64-v2":      (1, None),
}


def _cpu_flags():
    """The flag set of the running CPU (empty if /proc/cpuinfo is unreadable)."""
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("flags"):
                    return set(line.split(":", 1)[1].split())
    except OSError:
        pass
    return set()


def _select_library(names, flags):
    """Rank the available libfmod variants for a CPU with `flags`.

    `names` are library stems found on disk (libfmod_znver4, libfmod, ...).  A
    variant is usable when the CPU has the flag its -march implies; the generic
    'libfmod' from a plain `make` is always usable but always last, since its
    ARCH is whatever the builder chose.  Returns the usable stems, best first.
    Kept separate from the file system so it can be tested on any machine.
    """
    ranked = []
    for name in names:
        if name == "libfmod":
            ranked.append((-1, name))
            continue
        if not name.startswith("libfmod_"):
            continue
        rank, need = _ARCH_REQUIREMENTS.get(name[len("libfmod_"):], (0, None))
        if need is None or need in flags:
            ranked.append((rank, name))
    return [name for _, name in sorted(ranked, key=lambda t: -t[0])]


def _load_library():
    """Pick the fastest libfmod the current CPU can run.

    `make auto` builds one libfmod_<target>.so per distinct CPU type found in the
    cluster's node inventory (`libs/detect_arch.sh`), and `make dispatch` builds
    the fixed znver3/znver4 pair.  Selecting here, at import, is what lets ONE
    checkout run on every node yet use AVX-512 only where it is safe -- an
    AVX-512 binary started on a Zen 3 node dies with SIGILL.  A plain `make`
    leaves a single libfmod.so, which is used as the fallback.
    """
    stems = sorted(p.stem for p in _LIB_DIR.glob("libfmod*.so"))
    for name in _select_library(stems, _cpu_flags()):
        return np.ctypeslib.load_library(name, str(_LIB_DIR))
    # Nothing usable found by name; let load_library raise a clear error.
    return np.ctypeslib.load_library("libfmod", str(_LIB_DIR))


_lib = _load_library()


def i64(value):
    """byref-wrapped 64-bit integer scalar for ctypes Fortran calls."""
    return byref(c_int64(value))


def dbl(value):
    """byref-wrapped double scalar for ctypes Fortran calls."""
    return byref(c_double(value))
