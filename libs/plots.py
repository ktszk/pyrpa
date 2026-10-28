#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Figures of the main.py modes.  Everything that draws with matplotlib lives here, outside
libs.plibs (which holds the calculations only): band structure, 2D/3D Fermi surfaces, DOS,
spectral functions, BdG bands, dHvA frequencies, NMR, and the optical-response figures.
The functions take the arrays the plibs drivers return.
"""
import collections.abc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import libs.plibs as plibs


COLOR_GAP = 3  # color_option of the Eilenberger pairing-gap coloring (main.ColorMode.GAP)


def plot_dielectric(w: np.ndarray, eps: np.ndarray, fname: str = 'dielectric',
                    wp_scr: np.ndarray | None = None) -> None:
    """Publication-quality figure of the dielectric function eps(w).

    Journal conventions: single-column width, serif/mathtext labels, inward
    ticks on all four sides with minor ticks, thin recessive frame, and the two
    curves separated by BOTH colour and line style (Okabe-Ito CVD-safe pair +
    solid/dashed) so they stay distinct in colour-blind view and B/W print.

    @param     w: photon-energy grid [eV] (real, 1-D)
    @param   eps: complex dielectric function on w (eps1=Re, eps2=Im)
    @param fname: output basename; writes <fname>.pdf (vector) and <fname>.png
    @param wp_scr: screened plasma frequency [eV]; if given, marked with a guide
    Leaves the figure open so a later plt.show() displays it interactively.
    """
    blue, verm = '#0072B2', '#D55E00'          # Okabe-Ito colour-blind-safe pair
    rc = {'font.family': 'serif', 'mathtext.fontset': 'dejavuserif', 'font.size': 9,
          'axes.linewidth': 0.8, 'lines.linewidth': 1.3,
          'xtick.direction': 'in', 'ytick.direction': 'in',
          'xtick.top': True, 'ytick.right': True,
          'xtick.minor.visible': True, 'ytick.minor.visible': True,
          'legend.frameon': False, 'savefig.dpi': 600, 'figure.dpi': 600}
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(3.4, 2.6))   # PRB single-column width (~86 mm)
        ax.axhline(0.0, color='0.6', lw=0.6, zorder=0)
        if wp_scr is not None:
            # mark the screened plasma frequency (where eps1 crosses zero)
            ax.axvline(wp_scr, color='0.5', lw=0.8, ls=':', zorder=0)
            ax.annotate(r'$\hbar\omega_{\mathrm{p}}$', xy=(wp_scr, 0),
                        xytext=(3, 3), textcoords='offset points',
                        ha='left', va='bottom', fontsize=8, color='0.35')
        ax.plot(w, np.real(eps), color=blue, ls='-', label=r'$\varepsilon_1$ (real)')
        ax.plot(w, np.imag(eps), color=verm, ls='--', label=r'$\varepsilon_2$ (imag)')
        ax.set_xlabel(r'$\hbar\omega$ (eV)')
        ax.set_ylabel(r'Dielectric function $\varepsilon(\omega)$')
        ax.set_xlim(np.min(w), np.max(w))
        ax.legend(loc='best', handlelength=1.8)
        fig.tight_layout(pad=0.3)
        for ext in ('pdf', 'png'):
            fig.savefig(f'{fname}.{ext}', bbox_inches='tight')
    print(f"dielectric-function figure written to '{fname}.pdf' and '{fname}.png'", flush=True)


def plot_reflectivity(w: np.ndarray, refl: np.ndarray, fname: str = 'reflectivity', wp_scr: np.ndarray | None = None,
                      swatch: list | None = None) -> None:
    """Publication-quality figure of the normal-incidence reflectivity R(w).

    Same journal conventions as plot_dielectric (single-column width, serif
    mathtext, inward ticks, Okabe-Ito CVD-safe colours combined with distinct
    line styles). The three principal (xx, yy, zz) components are drawn; the
    yy/zz curves are omitted when they are numerically identical to xx (cubic
    or otherwise isotropic response), keeping the panel uncluttered.

    @param      w: photon-energy grid [eV] (real, 1-D)
    @param   refl: reflectivity array [Nw, 3] (xx, yy, zz), values in [0, 1]
    @param  fname: output basename; writes <fname>.pdf (vector) and <fname>.png
    @param wp_scr: screened plasma frequency [eV]; if given, marked with a guide
    @param swatch: optional [(label, rgb)] list; each is drawn as a colour patch
                   so the perceived metallic colour is shown next to R(w)
    Leaves the figure open so a later plt.show() displays it interactively.
    """
    blue, verm, green = '#0072B2', '#D55E00', '#009E73'   # Okabe-Ito CVD-safe triple
    sty = [(blue, '-', r'$R_{xx}$'), (verm, '--', r'$R_{yy}$'), (green, ':', r'$R_{zz}$')]
    rc = {'font.family': 'serif', 'mathtext.fontset': 'dejavuserif', 'font.size': 9,
          'axes.linewidth': 0.8, 'lines.linewidth': 1.3,
          'xtick.direction': 'in', 'ytick.direction': 'in',
          'xtick.top': True, 'ytick.right': True,
          'xtick.minor.visible': True, 'ytick.minor.visible': True,
          'legend.frameon': False, 'savefig.dpi': 600, 'figure.dpi': 600}
    with plt.rc_context(rc):
        fig, ax = plt.subplots(figsize=(3.4, 2.6))   # PRB single-column width (~86 mm)
        if wp_scr is not None:
            # the reflectivity edge sits at the screened plasma frequency
            ax.axvline(wp_scr, color='0.5', lw=0.8, ls=':', zorder=0)
            ax.annotate(r'$\hbar\omega_{\mathrm{p}}$', xy=(wp_scr, 1.0),
                        xytext=(3, -9), textcoords='offset points',
                        ha='left', va='top', fontsize=8, color='0.35')
        for i, (col, ls, lab) in enumerate(sty):
            # skip a component that duplicates xx (isotropic response)
            if i > 0 and np.allclose(refl[:, i], refl[:, 0], rtol=1e-10, atol=1e-12):
                continue
            ax.plot(w, refl[:, i], color=col, ls=ls, label=lab)
        ax.set_xlabel(r'$\hbar\omega$ (eV)')
        ax.set_ylabel(r'Reflectivity $R(\omega)$')
        ax.set_xlim(np.min(w), np.max(w))
        ax.set_ylim(0.0, 1.0)
        if swatch:
            # perceived colour patches (sRGB) as a small legend-like column in
            # the bottom-right corner: square patch + label, stacked upwards
            x0, y0, wsw, hsw, row = 0.560, 0.055, 0.042, 0.060, 0.075
            for j, (lab, rgb) in enumerate(reversed(swatch)):
                y = y0+j*row
                ax.add_patch(plt.Rectangle((x0, y), wsw, hsw, transform=ax.transAxes,
                                           facecolor=tuple(np.clip(rgb, 0, 1)),
                                           edgecolor='0.35', lw=0.5, zorder=5))
                ax.text(x0+wsw+0.018, y+0.5*hsw, lab, transform=ax.transAxes,
                        ha='left', va='center', fontsize=6, color='0.3', zorder=5)
        ax.legend(loc='best', handlelength=1.8)
        fig.tight_layout(pad=0.3)
        for ext in ('pdf', 'png'):
            fig.savefig(f'{fname}.{ext}', bbox_inches='tight')
    print(f"reflectivity figure written to '{fname}.pdf' and '{fname}.png'", flush=True)


def flatten_orbs(seq: list) -> list:
    """Flatten an olist (mix of ints and int-lists) into a flat list of orbital indices."""
    out = []
    for o in seq:
        out += o if isinstance(o, list) else [o]
    return out


def plot_band(eig: np.ndarray, spl: np.ndarray, xlabel: list, xticks: list, uni: np.ndarray, ol: np.ndarray,
              color: bool) -> None:
    """
    @brief Band structure along a symmetry line, optionally coloured by orbital weight
           (olist groups -> R, G, B).  All bands are drawn as one LineCollection: each segment
           k_i -> k_i+1 takes the colour of its starting point.
    @param    eig: band energies [Nband, Nk] (relative to mu)
    @param    spl: path coordinate [Nk]
    @param xlabel, xticks: labels and positions of the symmetry points
    @param    uni: eigenvectors [Nband, Norb, Nk]
    @param     ol: three orbital indices or index lists (R, G, B)
    @param  color: colour by orbital weight (False: black)
    """
    from matplotlib.collections import LineCollection

    def get_col(cl: np.ndarray, ol: np.ndarray) -> np.ndarray:
        col = (np.abs(cl[ol])**2 if isinstance(ol, int)
               else (np.abs(cl[ol])**2).sum(axis=0)).round(4)
        return col
    # Check if olist has required 3 elements
    if len(ol) < 3:
        print(f"Warning: plot_band olist requires 3 elements (current: {len(ol)}). Color display disabled", flush=True)
        color = False
    segs, cols = [], []
    for e, cl in zip(eig, uni):  # band loop
        pts = np.column_stack([spl, e])
        seg = np.stack([pts[:-1], pts[1:]], axis=1)          # [Nk-1, 2, 2]
        if color:
            # Normalize eigenvectors so that per-orbital weights (|u_n|^2) sum to 1
            norm = np.sqrt((abs(cl)**2).sum(axis=0))
            # Check for norm close to zero (use tolerance for floating-point comparison)
            if np.any(norm < 1e-14):
                print("Warning: plot_band found norm close to zero. Skipping", flush=True)
                continue
            cls = cl/norm
            # Map three orbital groups to RGB channels for color coding
            rgb = np.array([get_col(cls, ol[0]), get_col(cls, ol[1]), get_col(cls, ol[2])]).T
            cols.append(rgb[:-1])
        segs.append(seg)
    plt.figure()
    if segs:
        # round caps join the segments seamlessly, as the old one-plot-per-segment drawing did
        lc = LineCollection(np.concatenate(segs), linewidths=plt.rcParams['lines.linewidth'],
                            colors=np.concatenate(cols) if color else 'black', capstyle='round')
        plt.gca().add_collection(lc)
    for x in xticks[1:-1]:
        plt.axvline(x, ls='-', lw=0.25, color='black')
    plt.ylim(eig.min()*1.1, eig.max()*1.1)
    plt.xlim(0, spl.max())
    plt.axhline(0., ls='--', lw=0.25, color='black')
    plt.xticks(xticks, xlabel)
    plt.show()


def plot_FS(fscolors: list, klist: np.ndarray, color_option: int) -> None:
    fig = plt.figure()
    fig.add_subplot(111, aspect='equal')
    if color_option == 0:
        for kl in klist:
            for kk in kl:
                plt.plot(kk[:, 0], kk[:, 1], color='black', lw=2)
    else:
        if color_option == 2:
            v = []
            k = []
            for vl, kl in zip(fscolors, klist):
                for vv, kk in zip(vl, kl):
                    v.extend(vv)
                    k.extend(kk)
            v = np.array(v)
            k = np.array(k)
            vmax = v.max()
            vmin = v.min()
            # Check for vmax=vmin case
            if vmax == vmin:
                print(f"Warning: plot_FS has vmax=vmin={vmax}. Color map disabled", flush=True)
                vmax = vmin + 1.0  # Set default value
        for kl, fscol in zip(klist, fscolors):
            for kk, fcol in zip(kl, fscol):
                if color_option == 2:
                    clist = cm.jet((fcol-vmin)/(vmax-vmin))
                else:
                    clist = fcol
                for k1, k2, clst in zip(kk, kk[1:], clist):
                    plt.plot([k1[0], k2[0]], [k1[1], k2[1]], c=clst, lw=2)
        if color_option == 2:
            plt.scatter(k[:, 0], k[:, 1], s=0.1, c=v)
            plt.jet()
            plt.colorbar()
    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5, 0.5)
    plt.xticks([-0.5, 0, 0.5], [r'-$\pi$', '0', r'$\pi$'])
    plt.yticks([-0.5, 0, 0.5], [r'-$\pi$', '0', r'$\pi$'])
    plt.show()


def poly_vertices(polys: list) -> np.ndarray:
    """
    @brief Collect every vertex of a (possibly per-band nested) polygon collection.
           The polygons are ragged once the sheets are folded into the BZ -- clipping a
           triangle at a zone face leaves 3 to 6 vertices -- so they cannot simply be
           cast to one array.
    @param polys: polygons, or lists of polygons, with vertices in the last axis
    @return     : all vertices as one [N,3] array
    """
    out = []
    stack = list(polys)
    while stack:
        p = stack.pop()
        if isinstance(p, np.ndarray) and p.ndim >= 2 and p.shape[-1] == 3:
            out.append(p.reshape(-1, 3))
        else:
            stack.extend(list(p))
    return np.vstack(out) if out else np.zeros((0, 3))


def plot_3d_surf(fspolys: list, fscenters: np.ndarray, fscolors: list, surface_opt: int, bvec: np.ndarray) -> None:
    import mpl_toolkits.mplot3d.axes3d as axes3d
    import mpl_toolkits.mplot3d.art3d as art3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    clist = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'w']
    if surface_opt == 0:
        for j, polys in enumerate(fspolys):
            ax.add_collection3d(art3d.Poly3DCollection(polys, facecolor=clist[j%6]))
    else:
        if surface_opt == 1:
            tri = art3d.Poly3DCollection(fspolys, facecolors=fscolors, lw=0)
        else:
            if surface_opt == COLOR_GAP:  # diverging map, symmetric about 0 so the node sits on white
                cmap = cm.bwr
                vmax = np.abs(fscolors).max()
                vmin = -vmax
            else:
                cmap = cm.jet
                vmin, vmax = fscolors.min(), fscolors.max()
            # Check for vmax=vmin case
            if vmax == vmin:
                print(f"Warning: plot_3d_surf has vmax=vmin={vmax}. Using default colors", flush=True)
                nor_cols = np.ones_like(fscolors)*0.5
            else:
                nor_cols = (fscolors-vmin)/(vmax-vmin)
            tri = art3d.Poly3DCollection(fspolys, facecolors=cmap(nor_cols), lw=0)
            fs = ax.scatter(fscenters[:, 0], fscenters[:, 1], fscenters[:, 2],
                            c=fscolors, cmap=cmap, vmin=vmin, vmax=vmax, s=0.1)
            plt.colorbar(fs, format='%.2e')
        ax.add_collection3d(tri)
    ax.grid(False)
    bzv = BZedge(bvec, ax)
    # Sheets and BZ edges are both Cartesian k [A^-1] now, so one range fits both.
    # (kscale>1 samples beyond the 1st BZ, hence taking the max of the two extents.)
    verts = poly_vertices(fspolys)
    kmax = np.maximum(np.abs(bzv).max(axis=0), np.abs(verts).max(axis=0)) if len(verts) \
        else np.abs(bzv).max(axis=0)
    ax.set_xlim(-kmax[0], kmax[0])
    ax.set_ylim(-kmax[1], kmax[1])
    ax.set_zlim(-kmax[2], kmax[2])
    ax.set_box_aspect(tuple(kmax))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    plt.tight_layout()
    plt.show()


def set_init_3dfsplot(color_option: int, polys: list, centers: np.ndarray, blist: list, avec: np.ndarray,
                      rvec: np.ndarray, ham_r: np.ndarray, S_r: np.ndarray, olist: np.ndarray,
                      gap_sym: int | str | None = None, delta0: np.ndarray | None = None,
                      gap_orbital: np.ndarray | collections.abc.Callable | None = None,
                      gauge_kw: dict | None = None, ihbar: float = 1.0) -> tuple:
    if color_option == 0:
        fspolys = polys
        fscenters = []
        fscolors = []
        return fspolys, fscenters, fscolors
    else:
        if color_option == COLOR_GAP:  # the actual Eilenberger pairing gap, on the real 3D FS
            colors = plibs.gap_color_3d(centers, blist, rvec, ham_r, S_r, gap_sym, delta0, gap_orbital,
                                        **(gauge_kw or {}))
        else:
            colors = plibs.get_colors(centers, blist, ihbar*avec.T, rvec, ham_r, S_r, olist, color_option)
        fspolys = []
        fscenters = []
        fscolors = []
        for i in polys:
            fspolys.extend(i)
        for i in centers:
            fscenters.extend(i)
        for i in colors:
            fscolors.extend(i)
        # centers come back fractional; the sheets are Cartesian, so map the scatter
        # overlay with bvec to keep both in the same frame
        fscenters = np.array(fscenters).dot(plibs.get_bvec(avec))
        fscolors = np.array(fscolors)
        return fspolys, fscenters, fscolors


def plot_dos(wlist: np.ndarray, Dos: np.ndarray, olist: list, color_option: int) -> None:
    """
    @brief Plot the orbital-resolved DOS.
           color_option=1 draws the olist groups as red/green/blue lines over the total (black);
           color_option=0 stacks the olist groups as filled areas under the total.
    @param        wlist: energy mesh [Nw] (eV, relative to mu)
    @param          Dos: orbital-resolved DOS [Norb, Nw] (flibs.gen_dos)
    @param        olist: up to three orbital indices or index lists (see main.py olist)
    @param color_option: 0 (stacked areas) or 1 (orbital lines)
    """
    if color_option == 1:
        for ol, cl in zip(olist, ['red', 'green', 'blue']):
            if ol == []:
                continue
            if isinstance(ol, int):
                plt.plot(wlist, Dos[ol], color=cl)
            else:
                plt.plot(wlist, Dos[ol, :].sum(axis=0), color=cl)
        plt.plot(wlist, Dos.sum(axis=0), color='black')
    elif color_option == 0:
        clist = ['k', 'r', 'g', 'b', 'c', 'm', 'y']
        tmp_ol_b = flatten_orbs(olist)
        plt.fill_between(wlist, Dos[tmp_ol_b].sum(axis=0), Dos.sum(axis=0), color=clist[0])
        for i in range(len(olist)-1):
            tmp_ol_u = flatten_orbs(olist[i:])
            tmp_ol_b = flatten_orbs(olist[i+1:])
            plt.fill_between(wlist, Dos[tmp_ol_b].sum(axis=0), Dos[tmp_ol_u].sum(axis=0), color=clist[i+1])
        plt.fill_between(wlist, 0, Dos[tmp_ol_b].sum(axis=0), color=clist[len(olist)])
    plt.xlim(wlist[0], wlist[-1])
    plt.ylim(0, max(Dos.sum(axis=0))*1.2)
    plt.show()


def BZedge(bvec: np.ndarray, ax: object | None = None) -> np.ndarray:
    """
    @brief Draw first Brillouin zone boundary edges using Wigner-Seitz construction in reciprocal space.
           The edges are drawn in Cartesian k [Angstrom^-1], so whatever is overlaid on the same
           axes (Fermi sheets, k-points) must be Cartesian too.
    @param bvec: reciprocal lattice vectors as rows, shape (3,3), in Angstrom^-1
    @param   ax: matplotlib Axes3D; uses plt.gca() if None
    @return  the BZ corner points [Nv,3] in Cartesian k, for setting the axis ranges
    """
    if ax is None:
        ax = plt.gca()
    corners = []
    for face in plibs.bz_faces(bvec):
        poly = np.vstack([face, face[0]])
        ax.plot(poly[:, 0], poly[:, 1], poly[:, 2], color='black', lw=0.5)
        corners.append(face)
    return np.vstack(corners) if corners else np.zeros((0, 3))


def plot_path_map(spa_length: np.ndarray, xticks: list, xlabel: list, wlist: np.ndarray, Akw: np.ndarray) -> None:
    """
    @brief Colour map of a spectral function along a symmetry line (plibs.spectrum_on_path).
    @param spa_length: path coordinate of each k-point [Nk]
    @param  xticks, xlabel: positions and labels of the symmetry points
    @param      wlist: energies [Nw] (eV, relative to mu)
    @param        Akw: spectral function [Nk, Nw]
    """
    w, x = np.meshgrid(wlist, spa_length)
    plt.contourf(x, w, Akw, cmap=plt.hot())
    for x in xticks[1:-1]:
        plt.axvline(x, ls='-', lw=0.25, color='black')
    plt.xlim(0, spa_length.max())
    plt.axhline(0., ls='--', lw=0.25, color='black')
    plt.xticks(xticks, xlabel)
    plt.colorbar()
    plt.show()


def plot_sigma_iw(iwlist: np.ndarray, sigma_iw: np.ndarray) -> None:
    """
    @brief Quick look at one element of the self-energy on the Matsubara axis.
    @param   iwlist: Matsubara frequencies w_n [Nw] (eV)
    @param sigma_iw: Sigma_11(Gamma, iw_n) [Nw]
    """
    plt.plot(iwlist, sigma_iw.real, color='r', label='Re')
    plt.plot(iwlist, sigma_iw.imag, color='b', label='Im')
    plt.xlabel(r'$\omega_n$ (eV)')
    plt.ylabel(r'$\Sigma_{11}(\Gamma, i\omega_n)$')
    plt.legend()
    plt.show()


def plot_cpa_spectrum(spa_length: np.ndarray, xticks: list, xlabel: list, wlist: np.ndarray, Akw: np.ndarray,
                      fname: str = 'cpa_spectrum.png') -> None:
    """
    @brief CPA spectral function along the symmetry line (plibs.calc_cpa), saved to fname.
    @param spa_length: path coordinate of each k-point [Nk]
    @param  xticks, xlabel: positions and labels of the symmetry points
    @param      wlist: energies [Nw] (eV, relative to mu)
    @param        Akw: spectral function [Nk, Nw]
    @param      fname: output image
    """
    w, x = np.meshgrid(wlist, spa_length)
    plt.contourf(x, w, Akw, 100, cmap=plt.get_cmap('hot'))
    for xt in xticks[1:-1]:
        plt.axvline(xt, ls='-', lw=0.25, color='white')
    plt.xlim(0, spa_length.max())
    plt.axhline(0., ls='--', lw=0.25, color='white')
    plt.xticks(xticks, xlabel)
    plt.ylabel('Energy (eV)')
    plt.title('CPA spectrum')
    plt.colorbar(label=r'$A(\mathbf{k},\omega)$')
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.show()


def plot_susceptibility_path(sp: np.ndarray, w: np.ndarray, sus: np.ndarray, fname: str) -> None:
    """
    @brief Colour map |Im chi(q, w)| along a q-path (calc_path_spectrum, chis_spectrum_sc), saved to fname.
    @param sp, w: meshgrids of the path coordinate and the energy [Nq, Nw]
    @param   sus: susceptibility [Nq, Nw] complex
    @param fname: output image
    """
    plt.contourf(sp, w, abs(sus.imag), 100, cmap='hot')
    plt.colorbar()
    plt.savefig(fname=fname, dpi=300)
    plt.close()


def plot_qmap(qx: np.ndarray, qy: np.ndarray, z: np.ndarray, fname: str) -> None:
    """
    @brief Colour map of a real quantity on the (qx, qy) plane (chis_qmap / phi_qmap), saved to fname.
    @param qx, qy: meshgrids of the q-plane
    @param      z: values [Nqx, Nqy] (real; the caller takes |Re| or |Im|)
    @param  fname: output image
    """
    plt.contourf(qx, qy, z, 100, cmap='jet')
    plt.colorbar()
    plt.savefig(fname=fname, dpi=300)
    plt.close()


def plot_chis_qpoint(wlist: np.ndarray, chis: np.ndarray, chis_orb: np.ndarray, fname: str | None = None) -> None:
    """
    @brief Im chi_s(q, w) at one q-point: the total, then each orbital component (chis_q_point[_sc]).
    @param    wlist: energies [Nw]
    @param     chis: total susceptibility [Nw] complex
    @param chis_orb: orbital components [Nw, Ncomp] complex
    @param    fname: None: show both figures; else save fname.png and fname_orb.png
    """
    for y, tag in ((chis.imag[:, None], ''), (chis_orb.imag, '_orb')):
        plt.plot(wlist, y)
        if fname is None:
            plt.show()
        else:
            plt.savefig(f'{fname}{tag}.png', dpi=150)
        plt.close()


def plot_bdg_path(ky: np.ndarray, eig_bdg: np.ndarray, fname: str = 'BdG_band.png') -> None:
    """
    @brief BdG bands along (0,0,0) -> (0,0.5,0) (plibs.bdg_bands_path), saved to fname.
    @param      ky: path coordinate k_y [Nk]
    @param eig_bdg: BdG eigenvalues [Nk, 2*Norb]
    @param   fname: output image
    """
    for n in range(eig_bdg.shape[1]):
        plt.plot(ky, eig_bdg[:, n], 'b-', lw=0.5)
    plt.axhline(0, ls='--', lw=0.5, color='gray')
    plt.xlabel(r'$k_y \; (2\pi/a)$')
    plt.ylabel('Energy (eV)')
    plt.title('BdG bands  (0,0,0) → (0,0.5,0)')
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.clf()


def plot_dhva(all_results: dict, phi: float, fname: str = 'dhva_band.png') -> None:
    """
    @brief dHvA frequencies F(theta) per band, the extrema of one band joined into orbit
           branches by continuity of F (plibs.link_dhva_branches), saved to fname.
    @param all_results: band index -> [N, 2] array of (theta, F) from plibs.get_dhva_band
    @param         phi: azimuthal angle of the sweep [deg] (title)
    @param       fname: output image
    """
    _, ax = plt.subplots()
    prop_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for ci, band_idx in enumerate(sorted(all_results.keys())):
        d = all_results[band_idx]
        color = prop_cycle[ci % len(prop_cycle)]
        theta_F = {th: sorted(d[d[:, 0] == th, 1]) for th in np.unique(d[:, 0])}
        for i, br in enumerate(plibs.link_dhva_branches(theta_F)):
            pts = np.array(br)
            label = f'Band {band_idx+1}' if i == 0 else '_nolegend_'
            ax.plot(pts[:, 0], pts[:, 1], '-o', color=color, markersize=3, label=label)
    ax.set_xlabel('theta (deg)')
    ax.set_ylabel('F (T)')
    ax.set_title(f'dHvA frequency vs field angle  phi={phi:.1f} deg')
    ax.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.show()


def plot_tdf(e: np.ndarray, tdf: np.ndarray) -> None:
    """
    @brief Transport distribution function (xx component) of the Boltzmann mode with sw_tdf.
    @param   e: energies relative to mu [Nw] (eV)
    @param tdf: transport distribution function [Nw, 3, 3]
    """
    plt.plot(e, tdf[:, 0, 0])
    plt.show()


def plot_optical(res: dict) -> None:
    """
    @brief Figures of the linear-response optics (plibs.calc_conductivity_lrt): the dielectric
           function and the reflectivity with the perceived-colour swatches (both saved), and an
           interactive view of sigma_xx(w).
    @param res: dict returned by calc_conductivity_lrt (wlist, sigma, epsilon, refl, wp_scr,
                cols, ncols)
    """
    wlist, sigma, wp_scr = res['wlist'], res['sigma'], res['wp_scr']
    cols, ncols = res['cols'], res['ncols']
    # publication-quality dielectric-function figure (xx component), saved to file
    plot_dielectric(wlist[1:], res['epsilon'][1:, 0, 0], fname='dielectric', wp_scr=wp_scr)
    # publication-quality reflectivity figure (principal components)
    # one patch per distinct colour (a cubic metal collapses to a single swatch)
    sw = []
    for lab, c in cols:
        for e in sw:
            if e[0] == c['hex']:
                e[1].append(lab)
                break
        else:
            sw.append([c['hex'], [lab], c['rgb']])
    swatch = [(f'{",".join(labs)}  {hx}', rgb) for hx, labs, rgb in sw]
    # add the non-metal colours only when the sample actually lets light through
    # at this thickness; a metal is opaque, so they would all be black clutter
    if ncols[0][1]['Y'] > 0.05:
        swatch += [(f'{lab}  {c[key]}', c['rgb_neg'] if key == 'hex_neg' else c['rgb'])
                   for lab, c, key in ncols]
    plot_reflectivity(wlist[1:], res['refl'][1:], fname='reflectivity', wp_scr=wp_scr, swatch=swatch)
    # quick interactive view of the optical conductivity
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$\omega$ (eV)')
    ax.set_ylabel(r'$\sigma_{xx}(\omega)$ (S/m)')
    ax.plot(wlist, sigma[:, 0, 0].real, label='Re')
    ax.plot(wlist, sigma[:, 0, 0].imag, label='Im')
    ax.legend()
    plt.show()


def plot_nmr(res: dict, fname: str = 'nmr_sc.png') -> None:
    """
    @brief Two-panel summary plot of the NMR sweep.
    Left : K_s/K_n vs T/Tc (linear) -- Yosida curve for a singlet, finite residual
           for a triplet d||z.
    Right: (1/T1T)_s/(1/T1T)_n vs T/Tc (log-log) -- a Hebel-Slichter peak just below
           Tc for a full isotropic gap, its suppression for a sign-changing (s+-) gap,
           and the low-T power law (1/T1 ~ T^3 for line nodes, ~T^5 for point nodes,
           exponential for a full gap) read off as the slope.
    @param   res: dict returned by plibs.nmr_sweep
    @param fname: output png file name
    """
    t = res['t_red']
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))

    ax[0].plot(t, res['K_ratio'], 'o-', color='tab:blue', ms=4)
    ax[0].axhline(1.0, ls=':', lw=0.8, color='gray')
    ax[0].axvline(1.0, ls='--', lw=0.8, color='gray')
    ax[0].set_xlabel(r'$T/T_c$')
    ax[0].set_ylabel(r'$K_s/K_n$')
    ax[0].set_title('Knight shift')
    ax[0].set_xlim(0, max(1.05, t.max()))
    ax[0].set_ylim(bottom=0)

    ax[1].loglog(t, res['R_ratio'], 'o-', color='tab:red', ms=4)
    ax[1].axhline(1.0, ls=':', lw=0.8, color='gray')
    ax[1].axvline(1.0, ls='--', lw=0.8, color='gray')
    tt = t[t > 0]
    if len(tt):
        # The power law is quoted for 1/T1, but the axis shows 1/(T1 T): line nodes give
        # 1/T1 ~ T^3, hence a slope of 2 here.  Anchored at the highest temperature shown.
        ref = (tt / tt.max()) ** 2
        ax[1].loglog(tt, ref, 'k--', lw=0.8, label=r'$1/T_1\propto T^{3}$ (line nodes)')
        ax[1].legend(fontsize=8)
    ax[1].set_xlabel(r'$T/T_c$')
    ax[1].set_ylabel(r'$(1/T_1T)_s/(1/T_1T)_n$')
    ax[1].set_title(r'spin-lattice relaxation')

    fig.tight_layout()
    fig.savefig(fname, dpi=200)
    plt.close(fig)


def plot_band_mass(res: dict, xlabel: list, fname: str = 'electron_mass.png') -> None:
    """
    @brief Inverse band mass along the path, m_e/m_path, for every band (plibs.band_mass_path).
           The inverse is plotted because it stays finite where a band is flat along the path;
           its sign tells electron-like (>0) from hole-like (<0) curvature.  Saved to fname.
    @param    res: dict returned by plibs.band_mass_path
    @param xlabel: labels of the symmetry points
    @param  fname: output image
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        inv = 1.0/res['m_path']
    _, ax = plt.subplots()
    for n in range(inv.shape[1]):
        ax.plot(res['spa_length'], inv[:, n], lw=1.0, label=f'band {n+1}')
    for x in res['xticks'][1:-1]:
        ax.axvline(x, ls='-', lw=0.25, color='black')
    ax.axhline(0., ls='--', lw=0.25, color='black')
    ax.set_xlim(0, res['spa_length'].max())
    ax.set_xticks(res['xticks'], xlabel)
    ax.set_ylabel(r'$m_e/m^*$ (along the path)')
    if inv.shape[1] <= 12:
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.show()


def plot_impurity_spectrum(res: dict, xlabel: list, fname: str = 'imp_spectrum.png') -> None:
    """
    @brief Unfolded impurity spectral function A(k,w) along the symmetry line
           (plibs.calc_impurity_spectrum), saved to fname.
    @param    res: dict returned by plibs.calc_impurity_spectrum
    @param xlabel: labels of the symmetry points
    @param  fname: output image
    """
    w, x = np.meshgrid(res['wlist'], res['spa_length'])
    plt.contourf(x, w, res['Akw'], 100, cmap=plt.get_cmap('hot'))
    for xt in res['xticks'][1:-1]:
        plt.axvline(xt, ls='-', lw=0.25, color='white')
    plt.xlim(0, res['spa_length'].max())
    plt.axhline(0., ls='--', lw=0.25, color='white')
    plt.xticks(res['xticks'], xlabel)
    plt.ylabel('Energy (eV)')
    plt.title(f"impurity spectrum ({res['Nimp']}/{res['Nsite']} sites)")
    plt.colorbar(label=r'$A(\mathbf{k},\omega)$')
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.show()
