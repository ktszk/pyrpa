"""
Frequency grid of the convolution kernels: the IR (intermediate representation) basis of sparse-ir
as an alternative to the uniform Matsubara grid.

ir_grid builds what the Fortran module freq_grid needs (fg_set_ir); set_freq_grid switches it.
The work-array layout of a stored function: positive frequencies j = 1..Nw at slots 1..Nw, their
negatives at Nw+j (fermion) / Nw+j-1 (boson; slot 1 is v = 0, slot 2Nw unused).
"""
import numpy as np
import libs.flibs as flibs

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy and the sibling packages imported above.
__all__ = ['ir_grid', 'set_freq_grid', 'active_freq_grid', 'boson_frequencies', 'fermion_frequencies',
           'uniform_fermion_frequencies', 'to_uniform_fermion', 'ir_fit_fermion', 'ir_eval_fermion']

_active = None   # the grid set_freq_grid installed (None: uniform Matsubara)


def _bases(beta: float, wmax: float, eps: float) -> tuple:
    import sparse_ir
    from sparse_ir import sve
    import warnings
    with warnings.catch_warnings():   # sparse-ir 1.x: harmless np.divide(where=) warnings with recent numpy
        warnings.simplefilter('ignore', UserWarning)
        K = sparse_ir.LogisticKernel(beta*wmax)
        try:
            # 1.x: svd_strat='accurate', since its 'fast' interpolative SVD calls
            # scipy.linalg.interpolative.seed, which recent SciPy no longer has
            res = sve.compute(K, eps=eps, svd_strat='accurate')
        except TypeError:   # 2.x (C backend): no svd_strat
            res = sve.compute(K, eps=eps)
        return (sparse_ir.FiniteTempBasis('F', beta, wmax, eps=eps, sve_result=res),
                sparse_ir.FiniteTempBasis('B', beta, wmax, eps=eps, sve_result=res))


def ir_grid(temp: float, wmax: float, eps: float = 1.0e-7) -> dict:
    """
    @brief IR frequency grid for the Fortran convolution kernels (pass it to set_freq_grid).
    The sampling frequencies are sparse-ir's defaults (symmetric about 0; the bosonic set holds v = 0),
    the positive ones of the shorter statistics padded so that both store the same Nw.  The tau points
    are the union of the fermionic and bosonic ones (symmetric about beta/2, for the reflection
    tau -> beta - tau of the bubble).  Transform matrices in the work-array convention:
      forward  (Ntau, 2Nw) = U(tau) . fit / T,   backward (2Nw, Ntau) = T Uhat . fit_tau,
    so forward -> product -> backward is the raw sum sum_m a_m b_{n-m} of the uniform grid.
    @param temp: k_B T [eV];  wmax: cutoff of the spectra [eV] (Lambda = wmax/T);  eps: basis accuracy
    @return dict: Nw, Ntau, wF, wB (stored frequency values), fwdF, bwdF, fwdB, bwdB, rev (1-based
            tau -> beta - tau), static_wt (raw tau=0 midpoint weights, fermion), dens (physical tau=0
            weights), nF, nB (reduced frequencies of the slots), tau, fitF, and the bases bF, bB
    """
    beta = 1.0/temp
    bF, bB = _bases(beta, wmax, eps)
    nF = np.asarray(bF.default_matsubara_sampling_points())
    nB = np.asarray(bB.default_matsubara_sampling_points())
    pF, pB = np.sort(nF[nF > 0]), np.sort(nB[nB > 0])
    Nw = max(len(pF), len(pB)+1)
    while len(pF) < Nw:
        pF = np.append(pF, pF[-1]+2)
    while len(pB)+1 < Nw:
        pB = np.append(pB, pB[-1]+2)
    sF = np.concatenate([pF, -pF])                      # 2Nw fermionic slots
    sB = np.concatenate([[0], pB, -pB])                 # 2Nw-1 bosonic slots
    tau = np.unique(np.round(np.concatenate([bF.default_tau_sampling_points(),
                                             bB.default_tau_sampling_points()]), 12))
    rev = np.array([np.argmin(np.abs(tau-(beta-t))) for t in tau])
    if np.abs(tau[rev]-(beta-tau)).max() > 1.0e-9*beta:
        raise RuntimeError('ir_grid: the tau sampling points are not symmetric about beta/2')
    Ntau = len(tau)
    UhF, UtF = bF.uhat(sF), bF.u(tau)                   # (L, 2Nw), (L, Ntau)
    UhB, UtB = bB.uhat(sB), bB.u(tau)
    fitF, fitB = np.linalg.pinv(UhF.T), np.linalg.pinv(UhB.T)
    fwdF = UtF.T @ fitF/temp
    bwdF = temp*UhF.T @ np.linalg.pinv(UtF.T)
    fwdB = np.zeros((Ntau, 2*Nw), dtype=np.complex128)
    bwdB = np.zeros((2*Nw, Ntau), dtype=np.complex128)
    fwdB[:, :2*Nw-1] = UtB.T @ fitB/temp
    bwdB[:2*Nw-1, :] = temp*UhB.T @ np.linalg.pinv(UtB.T)
    mid = 0.5*(bF.u(np.array([0.0]))[:, 0]-bF.u(np.array([beta]))[:, 0])   # (G(0+) + G(0-))/2 per l
    dens = mid @ fitF
    return dict(Nw=Nw, Ntau=Ntau, wF=pF*np.pi*temp, wB=np.concatenate([[0], pB])*np.pi*temp,
                fwdF=np.asfortranarray(fwdF), bwdF=np.asfortranarray(bwdF),
                fwdB=np.asfortranarray(fwdB), bwdB=np.asfortranarray(bwdB),
                rev=(rev+1).astype(np.int64), static_wt=np.ascontiguousarray(dens/temp), dens=np.ascontiguousarray(dens),
                nF=sF, nB=sB, tau=tau, fitF=fitF, bF=bF, bB=bB, temp=temp)


def set_freq_grid(g: dict | None) -> None:
    """@brief Switch the Fortran kernels to the IR grid g (from ir_grid), or back to the uniform
    Matsubara grid with None."""
    global _active
    if g is None:
        flibs.fg_set_matsubara()
    else:
        flibs.fg_set_ir(g)
    _active = g


def active_freq_grid() -> dict | None:
    """@brief The IR grid installed by set_freq_grid, or None on the uniform Matsubara grid."""
    return _active


def _stored(key: str, temp: float, Nw: int) -> np.ndarray | None:
    if _active is None:
        return None
    if _active['Nw'] != Nw or abs(_active['temp']-temp) > 1.0e-12*temp:
        raise ValueError(f'the active IR grid has Nw = {_active["Nw"]}, T = {_active["temp"]:g} eV, not '
                         f'Nw = {Nw}, T = {temp:g} eV')
    return np.array(_active[key])


def boson_frequencies(temp: float, Nw: int) -> np.ndarray:
    """@brief The Nw stored bosonic frequencies v_j [eV] of the active grid (j = 0 is v = 0): 2 pi j T on the
    uniform Matsubara grid, the sampling frequencies on the IR grid.  For quantities that the caller builds
    as a function of v and hands to the kernels (e.g. the phonon vertex vph)."""
    w = _stored('wB', temp, Nw)
    return 2.0*np.pi*temp*np.arange(Nw) if w is None else w


def uniform_fermion_frequencies(temp: float, n: int) -> np.ndarray:
    """@brief The n lowest uniform fermionic Matsubara frequencies w_n = (2n+1) pi T [eV], whatever grid is active."""
    return (2*np.arange(n)+1)*np.pi*temp


def fermion_frequencies(temp: float, Nw: int) -> np.ndarray:
    """@brief The Nw stored fermionic frequencies w_n > 0 [eV] of the active grid: (2n+1) pi T on the uniform
    Matsubara grid, the sampling frequencies on the IR grid."""
    w = _stored('wF', temp, Nw)
    return uniform_fermion_frequencies(temp, Nw) if w is None else w


def ir_fit_fermion(g: dict, x_pos: np.ndarray, x_neg: np.ndarray | None = None) -> np.ndarray:
    """
    @brief IR coefficients x_l [..., L] of a fermionic function stored on the Nw positive sampling
    frequencies (last axis).  x_neg: its values at the negative ones; default conj(x_pos), i.e. a
    function with x(-iw) = x(iw)^* (a diagonal Green's-function element, Sigma_aa, ...).
    """
    x_pos = np.asarray(x_pos)
    x_neg = np.conj(x_pos) if x_neg is None else np.asarray(x_neg)
    return np.concatenate([x_pos, x_neg], axis=-1) @ g['fitF'].T


def ir_eval_fermion(g: dict, x_l: np.ndarray, n: np.ndarray) -> np.ndarray:
    """@brief Evaluate IR coefficients x_l [..., L] at fermionic Matsubara frequencies w = n pi T
    (n odd integers), e.g. the uniform w_j = (2j-1) pi T: n = 2j-1."""
    return x_l @ g['bF'].uhat(np.asarray(n))


def to_uniform_fermion(x: np.ndarray, temp: float, n_out: int) -> np.ndarray:
    """
    @brief A fermionic matrix function stored on the active grid, x [Norb, Norb, Nw, ...] (the layout of
    Sigma, G from flibs), evaluated at the n_out lowest uniform frequencies w_n = (2n+1) pi T.  On the
    uniform grid this is x[:, :, :n_out]; on the IR grid it goes through the basis, with the negative
    frequencies from x_ab(-iw) = x_ba(iw)^* (any Sigma or G, with or without SOC).
    The constant x(iw -> inf) (the Hartree-Fock part of Sigma) is not representable by the basis, so it is
    taken off first: the Hermitian part of x is even in 1/w, x(inf) + c2/w^2 + c4/w^4 + ..., and its values at
    the three highest sampling frequencies fix x(inf) up to O(1/w^6).
    @return [Norb, Norb, n_out, ...] complex128
    """
    x = np.asarray(x)
    if _active is None:
        if x.shape[2] < n_out:
            raise ValueError(f'to_uniform_fermion: only {x.shape[2]} frequencies stored, {n_out} requested')
        return x[:, :, :n_out]
    g = _active
    _stored('wF', temp, x.shape[2])                                     # Nw, T must be the grid's
    # Hermitian part at the three highest sampling frequencies = x(inf) + c2/w^2 + c4/w^4 (exact solve)
    herm = 0.5*(x[:, :, -3:]+np.conj(np.swapaxes(x[:, :, -3:], 0, 1)))
    iw2 = 1.0/g['wF'][-3:]**2
    c = np.linalg.solve(np.vander(iw2, 3, increasing=True), np.moveaxis(herm, 2, 0).reshape(3, -1))
    xinf = c[0].reshape(herm.shape[:2]+herm.shape[3:])[:, :, None]
    d = np.moveaxis(x-xinf, 2, -1)                                      # [..., Nw]
    dneg = np.conj(np.swapaxes(d, 0, 1))
    coef = ir_fit_fermion(g, d, dneg)
    out = ir_eval_fermion(g, coef, 2*np.arange(n_out)+1)                # [..., n_out]
    return np.moveaxis(out, -1, 2)+xinf
