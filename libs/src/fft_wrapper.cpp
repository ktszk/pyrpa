// fft_wrapper.cpp
// Drop-in FFTW3 API for Fortran via pocketfft (BSD-3-Clause).
// Provides dfftw_plan_dft_, dfftw_plan_many_dft_, dfftw_execute_,
// dfftw_execute_dft_ and dfftw_destroy_plan_ with Fortran name-mangled
// symbols (trailing underscore).
//
// Copyright (c) 2026 K. Suzuki
// SPDX-License-Identifier: MIT
// (see LICENSE; the bundled pocketfft_hdronly.h itself is BSD-3-Clause, see THIRD_PARTY_NOTICES)
// Speed-ups over a plain pocketfft call, both free of any change to the results' layout:
//  * the 1D plans (twiddle tables) of pocketfft are cached instead of rebuilt on every call;
//  * a transform called from serial code runs on several threads (pocketfft's own pool),
//    sized by OpenMP's setting; inside an OpenMP parallel region it stays single-threaded.
#ifndef POCKETFFT_CACHE_SIZE
#define POCKETFFT_CACHE_SIZE 16
#endif
#include "pocketfft_hdronly.h"
#include <complex>
#include <cstdint>
#include <cstddef>
#include <cstdlib>

using namespace pocketfft;

// OpenMP queries, bound WEAKLY to whichever OpenMP runtime the Fortran objects bring in
// (libgomp, LLVM libomp, libiomp5): linking one here could put two runtimes in the process.
// Without any runtime they resolve to null and the FFT runs single-threaded.
extern "C" {
int omp_get_max_threads(void) __attribute__((weak));
int omp_in_parallel(void) __attribute__((weak));
}

// Threads for one transform of ntot complex numbers.  PYRPA_FFT_THREADS overrides
// (1 = the old serial behaviour).  Otherwise one thread per 2^17 elements (2 MB), up to
// OpenMP's thread count: below that the pool's start-up costs more than it saves (measured
// on a 12-core Zen 4: 8x8x1x256 is 2x slower on 2 threads, 32x32x1x4096 is 4x faster on 24).
static size_t fft_threads(const shape_t &shape)
{
    if (const char *e = std::getenv("PYRPA_FFT_THREADS")) {
        const long v = std::strtol(e, nullptr, 10);
        if (v >= 1) return static_cast<size_t>(v);
    }
    if (omp_in_parallel && omp_in_parallel()) return 1;     // the caller already uses the cores
    if (!omp_get_max_threads) return 1;
    size_t ntot = 1;
    for (size_t n : shape) ntot *= n;
    const size_t tmax = static_cast<size_t>(omp_get_max_threads() > 1 ? omp_get_max_threads() : 1);
    const size_t twork = ntot >> 17;
    return twork < 1 ? 1 : (twork < tmax ? twork : tmax);
}

struct Plan {
    shape_t  shape;
    stride_t stride;       // input strides, in bytes
    stride_t stride_out;   // output strides, in bytes (equal to stride unless
                           // the advanced interface asked for a different layout)
    shape_t  axes;
    bool     forward;
    std::complex<double> *in;
    std::complex<double> *out;
};

extern "C" {

// Fortran signature:
//   integer(c_int64_t)          :: plan         (out, by reference)
//   integer(c_int32_t)          :: rank         (in,  by reference — literal 4)
//   integer(c_int32_t), dim(*) :: n            (in,  Nlist array)
//   complex(c_double), dim(*) :: in_arr, out_arr
//   integer(c_int32_t)          :: sign         (in,  -1=forward, +1=backward)
//   integer(c_int32_t)          :: flags        (in,  ignored)
void dfftw_plan_dft_(int64_t *plan_out,
                     const int32_t *rank,
                     const int32_t *n,
                     std::complex<double> *in_arr,
                     std::complex<double> *out_arr,
                     const int32_t *sign,
                     const int32_t * /*flags*/)
{
    auto *p = new Plan();
    int r = static_cast<int>(*rank);

    p->forward = (*sign == -1);   // FFTW_FORWARD=-1, FFTW_BACKWARD=+1
    p->in  = in_arr;
    p->out = out_arr;

    p->shape.resize(r);
    p->stride.resize(r);
    p->axes.resize(r);

    // Build column-major (Fortran) byte strides
    std::ptrdiff_t s = static_cast<std::ptrdiff_t>(sizeof(std::complex<double>));
    for (int i = 0; i < r; ++i) {
        p->shape[i] = static_cast<size_t>(n[i]);
        p->stride[i] = s;
        p->axes[i]   = static_cast<size_t>(i);
        s *= static_cast<std::ptrdiff_t>(n[i]);
    }
    p->stride_out = p->stride;

    *plan_out = static_cast<int64_t>(reinterpret_cast<uintptr_t>(p));
}

// FFTW advanced interface. Fortran signature:
//   dfftw_plan_many_dft(plan, rank, n, howmany,
//                       in,  inembed, istride, idist,
//                       out, onembed, ostride, odist,
//                       sign, flags)
// Modelled as a (rank+1)-dimensional array whose last dimension is the howmany
// batch: pocketfft loops over every axis that is absent from `axes`, which is
// exactly FFTW's batch semantics, so no explicit loop is needed here.
void dfftw_plan_many_dft_(int64_t *plan_out,
                          const int32_t *rank,
                          const int32_t *n,
                          const int32_t *howmany,
                          std::complex<double> *in_arr,
                          const int32_t *inembed,
                          const int32_t *istride,
                          const int32_t *idist,
                          std::complex<double> *out_arr,
                          const int32_t *onembed,
                          const int32_t *ostride,
                          const int32_t *odist,
                          const int32_t *sign,
                          const int32_t * /*flags*/)
{
    auto *p = new Plan();
    const int r = static_cast<int>(*rank);
    constexpr std::ptrdiff_t elem =
        static_cast<std::ptrdiff_t>(sizeof(std::complex<double>));

    p->forward = (*sign == -1);   // FFTW_FORWARD=-1, FFTW_BACKWARD=+1
    p->in  = in_arr;
    p->out = out_arr;

    // A null *embed means "contiguous, same as n" in FFTW.
    const int32_t *ne_in  = inembed ? inembed : n;
    const int32_t *ne_out = onembed ? onembed : n;

    p->shape.resize(r + 1);
    p->stride.resize(r + 1);
    p->stride_out.resize(r + 1);
    p->axes.resize(r);

    // Column-major (Fortran) layout: element (i_0..i_{r-1}) of transform j sits
    // at j*dist + sum_k i_k * stride * prod_{m<k} embed[m].
    std::ptrdiff_t si = static_cast<std::ptrdiff_t>(*istride) * elem;
    std::ptrdiff_t so = static_cast<std::ptrdiff_t>(*ostride) * elem;
    for (int i = 0; i < r; ++i) {
        p->shape[i]      = static_cast<size_t>(n[i]);
        p->stride[i]     = si;
        p->stride_out[i] = so;
        p->axes[i]       = static_cast<size_t>(i);
        si *= static_cast<std::ptrdiff_t>(ne_in[i]);
        so *= static_cast<std::ptrdiff_t>(ne_out[i]);
    }
    // Trailing batch axis, deliberately left out of `axes` so it is not transformed.
    p->shape[r]      = static_cast<size_t>(*howmany);
    p->stride[r]     = static_cast<std::ptrdiff_t>(*idist) * elem;
    p->stride_out[r] = static_cast<std::ptrdiff_t>(*odist) * elem;

    *plan_out = static_cast<int64_t>(reinterpret_cast<uintptr_t>(p));
}

void dfftw_execute_(const int64_t *plan)
{
    auto *p = reinterpret_cast<Plan *>(static_cast<uintptr_t>(*plan));
    // fct=1.0: Fortran side handles normalization (/product(Nlist))
    c2c(p->shape, p->stride, p->stride_out, p->axes, p->forward,
        p->in, p->out, 1.0, fft_threads(p->shape));
}

// FFTW threading calls, accepted for API compatibility with the MKL/FFTW backend: this
// wrapper picks the thread count of every execute by itself (fft_threads), so they do nothing.
void dfftw_init_threads_(int32_t *iret) { *iret = 1; }
void dfftw_plan_with_nthreads_(const int32_t * /*nthreads*/) {}

void dfftw_destroy_plan_(int64_t *plan)
{
    delete reinterpret_cast<Plan *>(static_cast<uintptr_t>(*plan));
    *plan = 0;
}

// Execute plan with different input/output arrays (same shape/type as planning arrays).
// Equivalent to FFTW's dfftw_execute_dft — safe to call from multiple threads
// as long as each thread uses distinct in_arr/out_arr buffers (each such call then
// runs single-threaded, see fft_threads).
void dfftw_execute_dft_(const int64_t *plan,
                        std::complex<double> *in_arr,
                        std::complex<double> *out_arr)
{
    const auto *p = reinterpret_cast<const Plan *>(static_cast<uintptr_t>(*plan));
    c2c(p->shape, p->stride, p->stride_out, p->axes, p->forward,
        in_arr, out_arr, 1.0, fft_threads(p->shape));
}

} // extern "C"
