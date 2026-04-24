// fft_wrapper.cpp
// Drop-in FFTW3 API for Fortran via pocketfft (BSD-3-Clause).
// Provides dfftw_plan_dft_, dfftw_execute_, dfftw_destroy_plan_
// with Fortran name-mangled symbols (trailing underscore).
//
// SPDX-License-Identifier: BSD-3-Clause
#include "pocketfft_hdronly.h"
#include <complex>
#include <cstdint>
#include <cstddef>

using namespace pocketfft;

struct Plan {
    shape_t  shape;
    stride_t stride;
    shape_t  axes;
    bool     forward;
    std::complex<double> *in;
    std::complex<double> *out;
};

extern "C" {

// Fortran signature:
//   integer(int64)          :: plan         (out, by reference)
//   integer(int32)          :: rank         (in,  by reference — literal 4)
//   integer(int32), dim(*) :: n            (in,  Nlist array)
//   complex(real64), dim(*) :: in_arr, out_arr
//   integer(int32)          :: sign         (in,  -1=forward, +1=backward)
//   integer(int32)          :: flags        (in,  ignored)
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

    *plan_out = static_cast<int64_t>(reinterpret_cast<uintptr_t>(p));
}

void dfftw_execute_(const int64_t *plan)
{
    auto *p = reinterpret_cast<Plan *>(static_cast<uintptr_t>(*plan));
    // fct=1.0: Fortran side handles normalization (/product(Nlist))
    c2c(p->shape, p->stride, p->stride, p->axes, p->forward,
        p->in, p->out, 1.0);
}

void dfftw_destroy_plan_(int64_t *plan)
{
    delete reinterpret_cast<Plan *>(static_cast<uintptr_t>(*plan));
    *plan = 0;
}

// Execute plan with different input/output arrays (same shape/type as planning arrays).
// Equivalent to FFTW's dfftw_execute_dft — safe to call from multiple threads
// as long as each thread uses distinct in_arr/out_arr buffers.
void dfftw_execute_dft_(const int64_t *plan,
                        std::complex<double> *in_arr,
                        std::complex<double> *out_arr)
{
    const auto *p = reinterpret_cast<const Plan *>(static_cast<uintptr_t>(*plan));
    c2c(p->shape, p->stride, p->stride, p->axes, p->forward,
        in_arr, out_arr, 1.0);
}

} // extern "C"
