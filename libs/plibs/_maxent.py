#!/usr/bin/env python
#-*- coding:utf-8 -*-
"""
Maximum-entropy (MaxEnt) analytic continuation of the trace of a Matsubara
Green's function to the real-frequency axis, as a drop-in, spurious-pole-free
alternative to flibs.pade_with_trace (Thiele continued fraction) for wide
frequency/k windows.

Requires the optional 'ana_cont' package (MIT license, pip install ana_cont).
"""
import numpy as np


def _ensure_np_trapz():
    """ana_cont (numpy<2.0 API) calls np.trapz, removed in numpy>=2.0 in favor
    of the identical np.trapezoid; restore the old name just before ana_cont
    needs it, rather than monkeypatching numpy at import of this module."""
    if not hasattr(np, 'trapz'):
        np.trapz = np.trapezoid


def maxent_trace_spectrum(Gk_mats: np.ndarray, wn: np.ndarray, wlist: np.ndarray, temp: float,
                          noise_level: float = 1.0e-4, alpha_determination: str = 'chi2kink',
                          verbose: bool = False, **kwargs) -> np.ndarray:
    """
    @fn maxent_trace_spectrum
    @brief MaxEnt analytic continuation of Tr G(k,iw_n) -> -Im Tr G(k,w+i0+), one
    k-point at a time. Replaces flibs.pade_with_trace when the Thiele continued
    fraction produces spurious poles over wide frequency/k windows: MaxEnt is a
    regularized (entropy-penalized) fit rather than an exact interpolant, so it
    does not chase numerical noise into extra real-axis structure.
    @param Gk_mats: Matsubara Green's function [Norb, Norb, Niw, Nk] complex128
                    (same layout as flibs.pade_with_trace's arrayin); only the
                    orbital-diagonal elements are used (their sum is Tr G, itself
                    a causal quantity with a positive spectral function)
    @param       wn: Matsubara frequencies [Niw] float64, e.g. (2n+1)*pi*temp
    @param    wlist: Real-frequency grid [Nw] float64 (no broadening term needed;
                     unlike Pade, MaxEnt does not require w+i*delta)
    @param     temp: Temperature in eV (beta = 1/temp)
    @param noise_level: Assumed relative noise of the input data (fraction of
                    max|Tr G|), used as the MaxEnt stdev; there is no actual
                    measurement error for FLEX/RPA output, so this stands in for
                    the numerical noise floor (self-consistency tolerance,
                    finite k/Matsubara mesh). Larger -> stronger regularization.
    @param alpha_determination: MaxEnt alpha-selection method (see ana_cont);
                    'chi2kink' is the library's recommended default.
    @param  verbose: forwarded to ana_cont (per-k, per-alpha-step logging)
    @return      Gk: -Im Tr G(k, w) [Nk, Nw] float64, same sign/scale convention
                     as -flibs.pade_with_trace(...).imag
    """
    _ensure_np_trapz()
    import ana_cont.continuation as cont

    Norb, _, Niw, Nk = Gk_mats.shape
    G_tr = np.einsum('iiwk->wk', Gk_mats)  # Tr G(iw_n, k); causal, positive spectral function

    beta = 1.0 / temp
    dw = wlist[-1] - wlist[0]
    model = np.full(wlist.shape, Norb / dw)  # flat default model, sum rule: int A dw = Norb
    stdev = np.full(Niw, noise_level * np.abs(G_tr).max())

    Gk = np.empty((Nk, len(wlist)), dtype=np.float64)
    for ik in range(Nk):
        problem = cont.AnalyticContinuationProblem(
            im_axis=wn, re_axis=wlist, im_data=G_tr[:, ik],
            kernel_mode='freq_fermionic', beta=beta)
        sol = problem.solve(method='maxent_svd', model=model, stdev=stdev,
                            alpha_determination=alpha_determination,
                            verbose=verbose, **kwargs)
        result = sol[0] if isinstance(sol, tuple) else sol
        Gk[ik, :] = np.pi * result.A_opt  # match -Im G scale of the Pade path
    return Gk
