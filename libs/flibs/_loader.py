import ctypes
import os
import pathlib
import numpy as np

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'i64', 'dbl', 'c_array', 'f_array'
]

# Resolve lib path relative to this file so import works from any current working directory.
_LIB_DIR = pathlib.Path(__file__).resolve().parents[1]


# -march target -> (rank, CPU flags it requires).  Higher rank = newer ISA, so the
# best library a CPU can run is the highest-ranked entry all of whose flags it has
# (the /proc/cpuinfo names).  `make auto` emits gcc's znver* names on AMD and the
# x86-64 levels on Intel (detect_arch.sh); the Intel CPU names are for an explicit
# `make ARCH=<name>`.  An unknown suffix is treated as rank 0 with no requirement,
# i.e. usable but least preferred.
_V3 = ("avx2", "fma", "bmi1", "bmi2", "f16c", "movbe")
_V4 = _V3 + ("avx512f", "avx512bw", "avx512cd", "avx512dq", "avx512vl")
_ARCH_REQUIREMENTS = {
    "znver5":         (5, ("avx512f",)),
    "znver4":         (4, ("avx512f",)),
    "x86-64-v4":      (4, _V4),
    "skylake-avx512": (4, _V4),
    "cascadelake":    (4, _V4 + ("avx512_vnni",)),
    "cooperlake":     (4, _V4 + ("avx512_vnni", "avx512_bf16")),
    "icelake-server": (4, _V4 + ("avx512_vnni", "avx512_vbmi2", "avx512_bitalg")),
    "icelake-client": (4, _V4 + ("avx512_vnni", "avx512_vbmi2", "avx512_bitalg")),
    "tigerlake":      (4, _V4 + ("avx512_vnni", "avx512_vbmi2", "avx512_vp2intersect")),
    "sapphirerapids": (4, _V4 + ("avx512_vnni", "avx512_vbmi2", "avx512_bf16", "amx_tile")),
    "emeraldrapids":  (4, _V4 + ("avx512_vnni", "avx512_vbmi2", "avx512_bf16", "amx_tile")),
    "graniterapids":  (4, _V4 + ("avx512_vnni", "avx512_vbmi2", "avx512_bf16", "amx_tile")),
    "znver3":         (2, ("avx2",)),
    "znver2":         (2, ("avx2",)),
    "znver1":         (2, ("avx2",)),
    "x86-64-v3":      (2, _V3),
    "haswell":        (2, _V3),
    "broadwell":      (2, _V3 + ("adx",)),
    "skylake":        (2, _V3 + ("adx",)),
    "alderlake":      (2, _V3 + ("adx", "avx_vnni")),
    "raptorlake":     (2, _V3 + ("adx", "avx_vnni")),
    "meteorlake":     (2, _V3 + ("adx", "avx_vnni")),
    "x86-64-v2":      (1, ("sse4_2", "popcnt")),
}


def _cpu_flags() -> set:
    """The flag set of the running CPU (empty if /proc/cpuinfo is unreadable)."""
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("flags"):
                    return set(line.split(":", 1)[1].split())
    except OSError:
        pass
    return set()


def _select_library(names: list, flags: set) -> list:
    """Rank the available libfmod variants for a CPU with `flags`.

    `names` are library stems found on disk (libfmod_znver4, libfmod, ...).  A
    variant is usable when the CPU has every flag its -march implies; the generic
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
        rank, need = _ARCH_REQUIREMENTS.get(name[len("libfmod_"):], (0, ()))
        if all(f in flags for f in need):
            ranked.append((rank, name))
    return [name for _, name in sorted(ranked, key=lambda t: -t[0])]


def _gfortran_runtime_env(path: pathlib.Path) -> None:
    """Environment a gfortran build needs, set before it is loaded unless the user chose otherwise.
    * OMP_WAIT_POLICY=PASSIVE: the build runs on libgomp, while numpy (MKL) keeps its own libiomp5
      pool in the same process; both spin-wait between parallel regions and fight for the cores,
      which slows a run down by ~100x on a fully used node.  libgomp reads it when it is loaded.
    * GFORTRAN_UNBUFFERED_PRECONNECTED=y: libgfortran buffers unit 6 when it is not a terminal and
      writes it only at exit, so the kernels' progress lines (the FLEX log the AiiDA steps read,
      'eliash=') would not reach a redirected stdout while the run goes on."""
    try:
        data = path.read_bytes()
    except OSError:
        return
    if b"libgomp.so" in data:
        os.environ.setdefault("OMP_WAIT_POLICY", "PASSIVE")
    if b"libgfortran.so" in data:
        os.environ.setdefault("GFORTRAN_UNBUFFERED_PRECONNECTED", "y")


def _load_library() -> np.ndarray:
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
        _gfortran_runtime_env(_LIB_DIR / f"{name}.so")
        return np.ctypeslib.load_library(name, str(_LIB_DIR))
    # Nothing usable found by name; let load_library raise a clear error.
    return np.ctypeslib.load_library("libfmod", str(_LIB_DIR))


_lib = _load_library()


def c_array(dtype: type) -> type:
    """
    @brief ctypes argtype for an array the Fortran side reads as one contiguous block in
    C order (the reversed-index view of a Fortran array, the convention of every wrapper).
    A plain ndpointer only checks the dtype and would hand Fortran the raw buffer of a
    strided view (a slice such as G[:, :, :Niw, :], or a transpose) without complaint, so
    the layout is enforced here: a non-contiguous argument raises TypeError at the call.
    @param dtype: numpy dtype of the array
    @return ndpointer class requiring C-contiguous memory
    """
    return np.ctypeslib.ndpointer(dtype=dtype, flags='C_CONTIGUOUS')


def f_array(dtype: type) -> type:
    """
    @brief ctypes argtype for an array passed in FORTRAN order (e.g. the orbital-pair list
    ol(Nchi,2), Umat(Norb,Norb) or mrot(3,3)): the Python index order equals the Fortran one.
    The argument is converted with np.asfortranarray at the call, so the result does not depend
    on the memory layout the caller happened to build (a C-ordered olist used to be read with
    its pairs scrambled).  An array that is already Fortran-contiguous is passed without a copy,
    so output arrays the wrappers allocate in Fortran order are still written in place.
    @param dtype: numpy dtype of the array
    @return ndpointer class taking any layout and handing Fortran-contiguous memory to the call
    """
    base = np.ctypeslib.ndpointer(dtype=dtype, flags='F_CONTIGUOUS')

    class _FortranOrder(base):
        @classmethod
        def from_param(cls, obj: object) -> object:
            # the returned ctypes view keeps the (possibly converted) array alive for the call
            return base.from_param(np.asfortranarray(obj, dtype=dtype))
    return _FortranOrder


def i64(value: int) -> object:
    """
    @brief byref-wrapped 64-bit integer scalar for ctypes Fortran calls.
    @param value: integer passed to Fortran by reference
    @return ctypes byref object pointing at a fresh c_int64
    """
    return ctypes.byref(ctypes.c_int64(value))


def dbl(value: float) -> object:
    """
    @brief byref-wrapped double scalar for ctypes Fortran calls.
    @param value: real number passed to Fortran by reference
    @return ctypes byref object pointing at a fresh c_double
    """
    return ctypes.byref(ctypes.c_double(value))
