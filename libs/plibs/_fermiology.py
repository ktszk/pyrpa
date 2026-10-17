#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Quantum-oscillation fermiology: cyclotron masses of extremal orbits and the dHvA
frequency F(theta) from the Onsager relation, on the Wannier-band Fermi surface, and the
band (effective) mass tensor along a symmetry line.
"""
import numpy as np
import libs.flibs as flibs
import libs.plibs._units as _units
import libs.plibs._bands as _bands
import libs.plibs._lattice as _lattice

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'get_mass', 'get_dhva_band', 'link_dhva_branches', 'band_mass_path', 'report_band_mass'
]


def get_mass(mesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, bvec: np.ndarray,
             mu: float, rotmat: np.ndarray, units: _units.Units, de: float = 3.e-6, meshkz: int = 20) -> None:
    """
    calculate cyclotron mass m*_c=hbar^2/2pi (dS/dE) where S is the area of Fermi surface cross section
    avec/bvec: primitive/reciprocal lattice vectors as rows, used for the slice area
    rotmat: rotation applied to the fractional R-vectors; its kz slices are the orbit planes
               Jacobian and to report the field direction rotmat actually describes
    de: energy step for numerical derivative (eV)
    meshkz: number of kz points for coarse scan of S(kz) in [0, 0.5] (reduced)
    """
    # Parameter validation
    if de <= 0:
        print("Error: Energy step (de) is non-positive", flush=True)
        return
    if meshkz <= 0:
        print("Error: kz mesh size (meshkz) is non-positive", flush=True)
        return

    # Area Jacobian of the kz-slice spanned by rotmat: |v1 x v2| with v_i = sum_k rotmat[i,k] b_k.
    # Equals 4*pi^2/(a1*a2) for an orthorhombic cell at rotmat=I, but follows the rotation
    # and the cell shape otherwise.
    ABZ = _bands.slice_area_factor(rotmat, bvec)
    bhat = _bands.field_direction(rotmat, avec)
    th_B = np.degrees(np.arccos(min(abs(bhat[2]), 1.)))
    ph_B = np.degrees(np.arctan2(bhat[1], bhat[0]))
    print(f"field direction implied by RotMat: B_hat = "
          f"({bhat[0]:.4f}, {bhat[1]:.4f}, {bhat[2]:.4f})  "
          f"theta={th_B:.2f} deg, phi={ph_B:.2f} deg", flush=True)
    print(f"slice area Jacobian = {ABZ:.6f} AA^-2", flush=True)
    _bands.open_contour_report()   # clear the counter for this run

    # --- Phase 1: coarse scan S(kz) ---
    print("Phase 1: scanning S(kz)...", flush=True)
    S_scan, eig_cache = _bands.scan_fs_area(mesh, rvec, ham_r, S_r, rotmat, mu, ABZ, meshkz)

    # --- Phase 2: refine extremal kz, compute m* ---
    print("Phase 2: computing m* at extremal orbits...", flush=True)
    print("="*50, flush=True)
    for band_idx in sorted(S_scan.keys()):
        data = np.array(S_scan[band_idx])
        kz_arr, S_arr = data[:, 0], data[:, 1]
        if len(kz_arr) < 3:
            continue
        cand_kz = _bands.find_extremal_kz(kz_arr, S_arr, band_idx, mesh, rvec, ham_r, S_r, rotmat, mu, ABZ)
        print(f"Band {band_idx+1}: extremal kz = {[f'{k:.4f}' for k in cand_kz]}", flush=True)

        results = []
        for kz_ext in cand_kz:
            eig = eig_cache.get(kz_ext)
            if eig is None:
                eig = _bands.get_eigs_2d(mesh, rvec, ham_r, S_r, rotmat, kz_ext)
            v2, blist = _bands.get_kf_points(eig, mesh, mu, kz_ext)
            S0 = _bands.get_band_area(v2, blist, band_idx, ABZ)
            if not S0:
                continue
            v2_p, blist_p = _bands.get_kf_points(eig, mesh, mu+de, kz_ext)
            v2_m, blist_m = _bands.get_kf_points(eig, mesh, mu-de, kz_ext)
            Sp = _bands.get_band_area(v2_p, blist_p, band_idx, ABZ)
            Sm = _bands.get_band_area(v2_m, blist_m, band_idx, ABZ)
            if Sp is None or Sm is None:
                continue
            # dS/dE via central finite difference; *1e20 converts AA^-2/eV -> m^-2/J
            dSdE_SI = (Sp-Sm)/(2.*de)*1.e20
            # Onsager: m*_c = (hbar^2 / 2pi) * dS/dE  [in units of electron mass m_e]
            mc = np.abs(dSdE_SI)*units.eC*units.hbar**2/(2*np.pi*units.emass)
            results.append((kz_ext, S0, mc))
            print(f"  kz={kz_ext:.4f}: S={S0:.4f} AA^-2, m*={mc:.4f} m_e", flush=True)

        if results:
            mc_vals = np.array([r[2] for r in results])
            kz_vals = np.array([r[0] for r in results])
            i_max, i_min = np.argmax(mc_vals), np.argmin(mc_vals)
            print(f"  >> Max m* = {mc_vals[i_max]:.4f} m_e  at kz = {kz_vals[i_max]:.4f}", flush=True)
            print(f"  >> Min m* = {mc_vals[i_min]:.4f} m_e  at kz = {kz_vals[i_min]:.4f}", flush=True)
    n_open = _bands.open_contour_report()
    if n_open:
        print(f"Warning: {n_open} orbit(s) discarded: the contour left the 2D slice box, so its "
              "area is not the orbit area. Increase the cell along the slice normal or use a "
              "field direction commensurate with the lattice.", flush=True)


def get_dhva_band(mesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, avec: np.ndarray, bvec: np.ndarray,
                  mu: float, theta_list: np.ndarray, units: _units.Units, phi: float = 0., meshkz: int = 20,
                  grid_mesh: int = 120) -> dict:
    """
    Calculate dHvA frequency F vs magnetic field polar angle theta.
    F = hbar/(2*pi*e) * A_ext  (Onsager relation)
    The Fermi surface is sliced by Cartesian planes perpendicular to B, over a window
    wide enough to close the orbit, so theta/phi may be any real angle: H(k) is exactly
    periodic in the fractional k it is fed, and the window may span several zones.
    avec/bvec : primitive/reciprocal lattice vectors as rows
    theta_list: array of polar angles from z-axis [deg]
    phi       : azimuthal angle [deg] (fixes the rotation plane, default 0 = xz-plane)
    (the figure is drawn by libs.plots.plot_dhva from the returned dict)
    meshkz    : number of slice offsets scanned per angle
    grid_mesh : diagonalise once on a grid_mesh^3 fractional grid and interpolate the
                slices from it (cubic). None diagonalises every slice point instead --
                exact but far slower; the grid reproduces it to ~1e-5 in the area.
    Returns   : dict  band_idx -> np.ndarray of shape (N, 2) columns=[theta, F[T]]
    """
    F_factor = units.hbar/(2.*np.pi)*1.e20   # AA^-2 -> T (Onsager): hbar[eV*s]/(2pi) * 1e20[AA^-2->m^-2]
    # keep the in-plane resolution the caller asked for: mesh points across a zone
    dk = float(np.linalg.norm(bvec, axis=1).min())/mesh
    egrid = None
    if grid_mesh:
        print(f"building the {grid_mesh}^3 band-energy grid...", flush=True)
        egrid = _bands.build_fs_energy_grid(rvec, ham_r, S_r, mu, mesh=grid_mesh)
        if egrid is None:
            print("Error: no band crosses the chemical potential", flush=True)
            return {}
        print(f"  bands kept: {[b+1 for b in egrid['bands']]}", flush=True)
    kf = _bands.fs_cartesian_points(rvec, ham_r, S_r, avec, mu, mesh=40, egrid=egrid)
    if len(kf) == 0:
        print("Error: no band crosses the chemical potential", flush=True)
        return {}

    all_results = {}  # band_idx -> [(theta, F), ...]
    has_axial = {}    # band_idx -> angles with a closed orbit encircling the field axis
    no_axial = {}     # band_idx -> angles without one
    seeds = []        # orbits found at earlier angles, followed into the next
    seed_life = 3     # angles a seed survives without being found again
    n_open_tot = 0

    for theta in theta_list:
        bhat = _bands.bfield_hat(theta, phi)
        radius, dmax = _bands.fs_extent_along(kf, bhat, avec)
        # commensurate B: A(d) repeats with period p and is symmetric about 0 and p/2,
        # so [0,p/2] is the fundamental domain and both ends are extremal planes.
        # Incommensurate B: no repeat, only interior dA/dd=0 orbits.
        period = _bands.slice_period(bvec, bhat)
        d_end = min(dmax, 0.5*period) if period is not None else dmax
        sym_end = d_end if (period is not None and 0.5*period <= dmax*(1.+1e-9)) else None
        L0, d_arr = 2.6*radius, np.linspace(0., d_end, meshkz, True)
        S_scan, orbit_cache, n_open, axial_bands = _bands.scan_plane_area(
            rvec, ham_r, S_r, avec, mu, bhat, d_arr, L0, dk, egrid=egrid)
        n_open_tot += n_open
        for b in S_scan:
            (has_axial if b in axial_bands else no_axial).setdefault(b, []).append(theta)
        # follow every orbit through the scan, take dA/dd=0 per branch, then collapse
        # the copies that differ by a whole reciprocal lattice vector
        branches = _bands.track_orbit_branches(orbit_cache, d_arr)
        ext = [e for br in branches for e in _bands.branch_extrema(br, d_end=sym_end)]
        # one step of the offset scan is the natural scale for telling a zone copy
        # (same orbit) from a neighbouring extremum of the same tube
        # Carry orbits found at earlier angles forward: an orbit already known is
        # followed straight to its extremum from its own centre, without having to be
        # rediscovered among hundreds of branches. This only ever adds. A seed that
        # fails to follow is kept for a few more angles rather than dropped at once --
        # one bad angle would otherwise break the chain for the whole rest of the sweep.
        step = d_end/max(meshkz-1, 1)
        alive = []
        for sd in seeds:
            f = _bands.follow_orbit(sd, rvec, ham_r, S_r, avec, mu, bhat, L0, dk, step, egrid=egrid)
            if f is None:
                if sd.get('_age', 0)+1 < seed_life:
                    sd = dict(sd); sd['_age'] = sd.get('_age', 0)+1
                    alive.append(sd)
                continue
            f['_age'] = 0
            ext.append(f)
            alive.append(f)
        ext = _bands.dedup_extremal_orbits(ext, bvec, tol_k=step, bhat=bhat)
        # an area that moves with the sampling window is not an orbit (see the docstring)
        ext, rejected = _bands.verify_orbits_window(ext, rvec, ham_r, S_r, avec, mu, bhat,
                                                   L0, dk, egrid=egrid)
        # seeds come from before the window check: continuity across angles is stronger
        # evidence than that check, and a rejected orbit must still be able to seed
        seeds = _bands.dedup_extremal_orbits(ext+alive, bvec, tol_k=step, bhat=bhat)
        for e in ext:
            all_results.setdefault(e['band'], []).append((theta, e['area']*F_factor))

        print(f"theta={theta:.1f} deg: FS bands={[b+1 for b in sorted(S_scan)]}, "
              f"window={L0:.2f} AA^-1, d in [0,{d_end:.3f}]"
              + ("" if period is None else f" (period {period:.3f})")+", "
              f"{len(branches)} branches -> {len(ext)} extremal orbits"
              + (f" ({len(rejected)} rejected: window dependent)" if rejected else "")
              + (f"  [no orbit encircles the axis: band(s) "
                 f"{[b+1 for b in sorted(set(S_scan)-axial_bands)]}]"
                 if set(S_scan)-axial_bands else ""), flush=True)

    # Only worth saying for a band that HAS a central orbit somewhere in the sweep:
    # a sheet that never encircles the axis (electron pockets off at X, say) is not an
    # open orbit, it simply lives elsewhere.
    for b in sorted(set(has_axial) & set(no_axial)):
        print(f"Band {b+1}: no closed orbit encircles the field axis for theta = "
              f"{', '.join(f'{t:g}' for t in no_axial[b])} deg, so the central orbit of "
              "this sheet gives no dHvA oscillation there (open orbits show up in "
              "magnetoresistance instead).", flush=True)
    if n_open_tot:
        print(f"Note: {n_open_tot} contour piece(s) never closed in the sampling window; "
              "open orbits and window-cut copies together.", flush=True)

    # convert to arrays
    all_results = {k: np.array(v) for k, v in all_results.items()}
    if not all_results:
        print("Error: no usable extremal orbit was found at any angle", flush=True)
    return all_results


def link_dhva_branches(theta_F: dict, tol: float = 0.35) -> list:
    """
    @brief Connect the extremal frequencies of one band into orbit branches across the field
           angle by continuity of F (a rank sort would swap branches wherever two of them cross).
    @param theta_F: theta -> sorted list of the frequencies F found at that angle
    @param     tol: largest relative jump |F-F'|/max(F,F') still linked into one branch
    @return list of branches, each a list of (theta, F)
    """
    out = []
    active = []                      # (branch index, last F)
    for th in sorted(theta_F):
        fs = theta_F[th]
        pairs = sorted(((abs(f-fp)/max(f, fp, 1e-30), ia, ifq)
                        for ia, (_, fp) in enumerate(active)
                        for ifq, f in enumerate(fs)), key=lambda t: t[0])
        used_a, used_f, new_active = set(), set(), []
        for cost, ia, ifq in pairs:
            if cost > tol or ia in used_a or ifq in used_f:
                continue
            used_a.add(ia); used_f.add(ifq)
            bi = active[ia][0]
            out[bi].append((th, fs[ifq]))
            new_active.append((bi, fs[ifq]))
        for ifq, f in enumerate(fs):
            if ifq not in used_f:
                out.append([(th, f)])
                new_active.append((len(out)-1, f))
        active = new_active
    return out


def band_mass_path(k_sets: list, kmesh: int, rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray,
                   avec: np.ndarray, units: _units.Units, degen_tol: float = 1.0e-4) -> dict:
    """
    @brief Band (effective) mass tensor along a symmetry line, from the analytic second
           derivative of the Wannier bands (intra- plus inter-band term, flibs.get_mass).
           The INVERSE mass tensor is the primary quantity: it stays finite for flat or
           strictly 2D directions (1/m = 0 there), where the mass itself diverges.
             inv_mass[k,n,a,b] = m_e (1/hbar^2) d^2 e_n/dk_a dk_b   (in units of 1/m_e)
           Three masses are derived from it, all in units of m_e (inf where 1/m = 0):
             m_path : along the path direction t at that k, 1/(t.M^-1.t) -- the curvature of
                      the plotted band; at a symmetry point t is the outgoing segment
             m_axes : along x, y, z, 1/(M^-1)_aa
             m_princ: principal masses, 1/eigenvalues of M^-1, sorted by eigenvalue
           A band degenerate with another one within degen_tol [eV] has no well-defined mass
           (the band index is ill-defined there): its entries are set to nan.
    @param k_sets, kmesh: symmetry points and points per segment (as for the band plot)
    @param rvec, ham_r, S_r: Wannier R-vectors, Hamiltonian, overlap ([] if orthogonal)
    @param   avec: primitive translation vectors as rows [Angstrom]
    @param  units: physical constants (unit_constants)
    @param degen_tol: degeneracy window [eV]
    @return dict: klist, spa_length, xticks, eig [Nk,Nb], inv_mass [Nk,Nb,3,3], m_path [Nk,Nb],
                  m_axes [Nk,Nb,3], m_princ [Nk,Nb,3]
    """
    bvec = _lattice.get_bvec(avec)
    klist, spa_length, xticks = _lattice.mk_klist(k_sets, kmesh, bvec)
    eig, uni = _bands.get_eigs(klist, ham_r, S_r, rvec)
    # get_mass with mrot = avec.T*ihbar gives (1/hbar^2) d^2e/dk^2 in eV m^2/(eV s)^2; dividing
    # by e converts eV -> J, i.e. to 1/kg, and multiplying by m_e to units of 1/m_e
    imass = flibs.get_mass(klist, ham_r, rvec, avec.T*units.ihbar, uni, True, S_r, eig)
    inv_mass = imass*units.emass/units.eC
    # path direction in Cartesian k: outgoing segment (incoming one at the last point)
    kc = np.asarray(klist).dot(bvec)
    dk = np.diff(kc, axis=0)
    dk = np.vstack([dk, dk[-1:]])
    nrm = np.linalg.norm(dk, axis=1)
    t = np.divide(dk, nrm[:, None], out=np.zeros_like(dk), where=nrm[:, None] > 0)

    def inv(x: np.ndarray) -> np.ndarray:
        with np.errstate(divide='ignore'):
            return np.where(np.abs(x) > 1.0e-12, 1.0/np.where(x == 0, 1.0, x), np.inf)
    m_path = inv(np.einsum('ka,knab,kb->kn', t, inv_mass, t))
    m_axes = inv(np.einsum('knaa->kna', inv_mass))
    m_princ = inv(np.linalg.eigvalsh(0.5*(inv_mass+np.swapaxes(inv_mass, 2, 3))))
    # degenerate bands: mask
    d = np.diff(eig, axis=1)
    deg = np.zeros(eig.shape, dtype=bool)
    deg[:, 1:] |= d < degen_tol
    deg[:, :-1] |= d < degen_tol
    for x in (m_path, m_axes, m_princ):
        x[deg] = np.nan
    return {'klist': klist, 'spa_length': spa_length, 'xticks': xticks, 'eig': eig, 'inv_mass': inv_mass,
            'm_path': m_path, 'm_axes': m_axes, 'm_princ': m_princ}


def report_band_mass(res: dict, xlabel: list, mu: float, fname: str = 'electron_mass.dat') -> None:
    """
    @brief Print the band masses at the symmetry points and write the whole path to fname.
    @param    res: dict returned by band_mass_path
    @param xlabel: labels of the symmetry points
    @param     mu: chemical potential [eV] (energies are reported relative to it)
    @param  fname: output file; one block per band: path coordinate, e-mu, m_path, m_x, m_y,
                   m_z, three principal masses (m_e; inf = flat direction, nan = degenerate)
    """
    eig = res['eig']-mu
    Nk = len(eig)
    # symmetry points sit at the start of every segment (kmesh points apart) and at the end
    nseg = len(xlabel)-1
    ipts = [round(i*(Nk-1)/nseg) for i in range(nseg+1)]
    print('band masses at the symmetry points (m_e; inf: flat direction, nan: degenerate band)', flush=True)
    print('  point  band    e-mu[eV]     m_path        m_x        m_y        m_z', flush=True)
    for lab, i in zip(xlabel, ipts):
        name = lab.replace('$', '').replace('\\', '')
        for n in range(eig.shape[1]):
            mx, my, mz = res['m_axes'][i, n]
            print(f'  {name:>5s} {n+1:5d} {eig[i, n]:+11.4f} {res["m_path"][i, n]:+10.4f} '
                  f'{mx:+10.4f} {my:+10.4f} {mz:+10.4f}', flush=True)
    try:
        with open(fname, 'w') as f:
            f.write('# path  e-mu[eV]  m_path  m_x  m_y  m_z  m_princ1  m_princ2  m_princ3  [m_e]\n')
            for n in range(eig.shape[1]):
                f.write(f'# band {n+1}\n')
                for i in range(Nk):
                    vals = [res['m_path'][i, n], *res['m_axes'][i, n], *res['m_princ'][i, n]]
                    f.write(f'{res["spa_length"][i]:10.6f} {eig[i, n]:+12.6f} '
                            + ' '.join(f'{v:+13.6e}' for v in vals) + '\n')
                f.write('\n')
        print(f"band masses along the path written to '{fname}'", flush=True)
    except OSError as e:
        print(f"Error: Failed to write '{fname}': {e}", flush=True)
