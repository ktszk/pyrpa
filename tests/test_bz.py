#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Tests for the 1st Brillouin zone geometry used by the 3D Fermi-surface plot
(libs/plibs/_lattice.py: bz_faces, bz_face_planes, fold_to_bz, and the folding
step inside libs/plibs/_bands.py: gen_3d_surf_points).

Marching cubes runs on the parallelepiped k_frac in [-1/2,1/2]^3, i.e. on the
primitive reciprocal CELL, while BZedge draws the Wigner-Seitz cell.  The two
coincide only for a primitive lattice; on a body/face/base-centred or otherwise
skewed cell the sheet sticks out of the drawn zone unless it is folded back, so
what has to be pinned is:

  * the zone really is the Wigner-Seitz cell -- 6 faces for simple cubic, 14
    (truncated octahedron) for fcc, 12 (rhombic dodecahedron) for bcc, 8 for
    hexagonal -- with no redundant tangent planes among them
  * folding moves nothing: the translates of the parallelepiped tile k-space, so
    the folded polygons cover the same total area, with nothing duplicated
  * after folding, no vertex of any sheet (nor any polygon centre fed back to
    get_colors) lies outside the zone
  * a primitive lattice is left strictly untouched

Runs standalone (no pytest needed):  python tests/test_bz.py
Also works under pytest if installed:  pytest tests/test_bz.py
Requires the Fortran library libfmod.so (cd libs && make).
"""
import sys
import os
import itertools

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

import libs.plibs as P

# primitive translation vectors of the lattices under test [Angstrom]
AVEC = {
    'cubic': 3.0*np.eye(3),
    'fcc': 3.6*np.array([[0., .5, .5], [.5, 0., .5], [.5, .5, 0.]]),
    'bcc': 3.0*np.array([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]),
    'hexagonal': np.array([[3., 0., 0.], [-1.5, 1.5*np.sqrt(3.), 0.], [0., 0., 10.]]),
}
NFACE = {'cubic': 6, 'fcc': 14, 'bcc': 12, 'hexagonal': 8}


def _poly_area(polys):
    """Total area of a list of convex polygons, summed over a triangle fan."""
    tot = 0.
    for p in polys:
        p = np.asarray(p, dtype=float)
        tot += .5*sum(np.linalg.norm(np.cross(b-p[0], a-p[0])) for a, b in zip(p[1:-1], p[2:]))
    return tot


def _nn_model(avec):
    """Single-orbital nearest-neighbour model on the given lattice."""
    rv = np.array([r for r in itertools.product((-1, 0, 1), repeat=3) if r != (0, 0, 0)], dtype=float)
    dist = np.linalg.norm(rv.dot(avec), axis=1)
    rvec = rv[np.isclose(dist, dist.min())]
    return rvec, np.full((len(rvec), 1, 1), -0.5+0j)


def _outside(points, bvec):
    """Largest signed distance of any point beyond a zone face (<=0 means inside)."""
    nvec, dist = P.bz_face_planes(bvec)
    return float((np.asarray(points, dtype=float).reshape(-1, 3).dot(nvec.T)-dist).max())


def test_zone_is_the_wigner_seitz_cell():
    """Face count and face distance of the 1st BZ, and no redundant planes: a bisector
    that only touches the cell along an edge (G=(1,1,0) for a cubic lattice) must not be
    counted as a face, or the clipper wastes passes on planes that cut nothing."""
    for name, avec in AVEC.items():
        bvec = P.get_bvec(avec)
        nvec, dist = P.bz_face_planes(bvec)
        assert len(nvec) == NFACE[name], f'{name}: {len(nvec)} faces, expected {NFACE[name]}'
        assert np.allclose(np.linalg.norm(nvec, axis=1), 1.), f'{name}: normals are not unit'
        assert len(P.bz_faces(bvec)) == NFACE[name], f'{name}: bz_faces disagrees with bz_face_planes'
        # every face polygon sits ON its own plane and inside all the others
        for face in P.bz_faces(bvec):
            s = face.dot(nvec.T)-dist
            assert np.isclose(s.max(), 0., atol=1e-9*dist.max()), f'{name}: a face lies outside the zone'
            assert np.isclose(np.abs(s).min(), 0., atol=1e-9*dist.max()), f'{name}: a face is on no plane'


def test_parallelepiped_sticks_out_of_the_zone():
    """The premise of the fold: the sampled cell is NOT the zone for a centred lattice."""
    corners = np.array(list(itertools.product((-.5, .5), repeat=3)))
    for name, avec in AVEC.items():
        bvec = P.get_bvec(avec)
        out = _outside(corners.dot(bvec), bvec)
        if name == 'cubic':
            assert out <= 1e-12, 'the cubic cell should already be the zone'
        else:
            assert out > 1e-3, f'{name}: the sampled cell was expected to overflow the zone'


def test_folding_conserves_area_and_lands_inside():
    """Folding is a bijection of the sampled cell onto the zone: no area is lost,
    duplicated, or left outside."""
    rng = np.random.default_rng(0)
    tris = rng.uniform(-.5, .5, size=(300, 3, 3))
    for name, avec in AVEC.items():
        bvec = P.get_bvec(avec)
        cart = tris.dot(bvec)
        folded = P.fold_to_bz(cart, bvec)
        assert _outside(np.vstack([np.asarray(p) for p in folded]), bvec) <= 1e-9, \
            f'{name}: a folded polygon is still outside the zone'
        assert np.isclose(_poly_area(folded), _poly_area(cart), rtol=1e-10), \
            f'{name}: folding changed the total area'


def test_primitive_lattice_is_untouched():
    """Nothing to fold on a simple cubic cell: the triangles must come back as they were."""
    rng = np.random.default_rng(1)
    bvec = P.get_bvec(AVEC['cubic'])
    cart = rng.uniform(-.5, .5, size=(50, 3, 3)).dot(bvec)
    folded = P.fold_to_bz(cart, bvec)
    assert len(folded) == len(cart)
    assert np.array_equal(np.array(folded), cart)


def test_fermi_sheet_is_folded_into_the_zone():
    """The reported bug: with a centred lattice the 3D Fermi surface left the BZ frame.
    The sheet must end up inside the zone, with its area (and the polygon centres that
    carry the colours) unchanged."""
    for name in ('fcc', 'bcc'):
        avec = AVEC[name]
        bvec = P.get_bvec(avec)
        rvec, ham_r = _nn_model(avec)
        S_r = np.array([])
        raw, _, braw = P.gen_3d_surf_points(16, rvec, ham_r, S_r, 1.0, 1.0, bvec, sw_fold=False)
        fold, cent, bfold = P.gen_3d_surf_points(16, rvec, ham_r, S_r, 1.0, 1.0, bvec)
        assert braw == bfold and len(raw) > 0, f'{name}: no Fermi sheet to test'
        assert _outside([v for p in raw for v in np.asarray(p).reshape(-1, 3)], bvec) > 1e-3, \
            f'{name}: the unfolded sheet was expected to leave the zone'
        for sheet, centre in zip(fold, cent):
            verts = np.vstack([np.asarray(p) for p in sheet])
            assert _outside(verts, bvec) <= 1e-9, f'{name}: the folded sheet leaves the zone'
            # the centres are fed back to get_eigs/get_colors as FRACTIONAL k
            assert len(centre) == len(sheet), f'{name}: one centre per polygon'
            assert _outside(np.asarray(centre).dot(bvec), bvec) <= 1e-9, \
                f'{name}: a polygon centre lies outside the zone'
        assert np.isclose(sum(_poly_area(s) for s in fold),
                          sum(_poly_area(np.asarray(s)) for s in raw), rtol=1e-8), \
            f'{name}: folding changed the Fermi-surface area'


def test_kscale_disables_the_fold():
    """kscale != 1 is a deliberate request for another k-window (the manual sells it as a
    display scale), so the sheet must be left where the sampling put it."""
    avec = AVEC['fcc']
    bvec = P.get_bvec(avec)
    rvec, ham_r = _nn_model(avec)
    sheets, _, _ = P.gen_3d_surf_points(12, rvec, ham_r, np.array([]), 1.0, [1., 1., .5], bvec)
    assert isinstance(sheets[0], np.ndarray) and sheets[0].shape[1:] == (3, 3), \
        'the sheet was folded although kscale != 1'


# --------------------------------------------------------------------------- #
#  standalone runner (no pytest required)
# --------------------------------------------------------------------------- #
if __name__ == '__main__':
    import _tools
    sys.exit(_tools.run_standalone(globals()))
