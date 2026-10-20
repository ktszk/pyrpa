"""
Python model calculation library (calculations only; the figures live in libs.plots).

Functions are organized in submodules by physical topic:
  _symmetry   : Point group operations, symmetry detection, irreducible k-points
  _lattice    : Lattice vectors, hopping import, k/r-mesh generation, BZ visualization
  _bands      : Eigenvalues, chemical potential, Fermi surface
  _response   : Spin susceptibility, pairing susceptibility, conductivity, gap symmetries
  _optical    : Optical constants (N, R, alpha) and the non-metal colour spectra (T, Kubelka-Munk)
  _color      : Colorimetry, perceived colour (sRGB) from a reflectivity/transmittance spectrum
  _nmr        : NMR observables in the SC state (Knight shift, 1/T1T) from chi_s
  _wannier_io : EPA output reading, U/J parameter sets, Wannier-R space file I/O
  _maxent     : MaxEnt analytic continuation (optional 'ana_cont' dependency)
  _eliash_nosym: Linearized Eliashberg without TRS/inversion (Zeeman field, FFLO pair momentum Q)
  _units      : Physical constants (SI/eV, or all 1 in the dimensionless test mode)
  _conductivity: Boltzmann and linear-response transport drivers (tau models, Hall, optics)
  _fermiology : Cyclotron mass, dHvA frequencies of extremal orbits, band masses on a path
  _impurity   : Impurity supercell, exact diagonalization and the unfolded spectrum
  _phonon_pairing: Phonon-mediated pairing interaction from EPA data for the Eliashberg equation
"""
from ._symmetry import *
from ._lattice import *
from ._bands import *
from ._response import *
from ._optical import *
from ._color import *
from ._nmr import *
from ._wannier_io import *
from ._maxent import *
from ._calc import *
from ._eliash_nosym import *
from ._eilenberger import *          # homogeneous solver + model Fermi surface (build_model_fs, ...)
from ._eilenberger_surface import *  # specular surface (gap profile, LDOS)
from ._eilenberger_vortex import *   # isolated vortex + periodic vortex lattice
from ._eilenberger_spin import *     # 2x2 spin-matrix (Pauli, d-vector textures)
from ._units import *
from ._conductivity import *
from ._fermiology import *
from ._impurity import *
from ._phonon_pairing import *
