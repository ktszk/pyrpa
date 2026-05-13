"""
Python model calculation library.

Functions are organized in submodules by physical topic:
  _symmetry   : Point group operations, symmetry detection, irreducible k-points
  _lattice    : Lattice vectors, hopping import, k/r-mesh generation, BZ visualization
  _bands      : Eigenvalues, chemical potential, Fermi surface
  _response   : Spin susceptibility, pairing susceptibility, conductivity, gap symmetries
  _wannier_io : EPA output reading, Wannier-R space file I/O
"""
from ._symmetry   import *
from ._lattice    import *
from ._bands      import *
from ._response   import *
from ._wannier_io import *
from ._calc     import *
