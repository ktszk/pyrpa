# Stored U/J parameter sets

Orbital-dependent Coulomb (`U`) and Hund (`J`) matrices, one TOML file per material, loaded
with `UJ_material` (see the manual, `UJ_material` / `UJ_dir`).

Each file carries its own reference (`source`, `title`, `doi`, `arxiv`), which `read_UJ`
prints when the set is loaded; please cite that paper when you use a set.

```toml
# comments are allowed
material = "FeSe"
source   = "authors, journal volume, page (year)"
title    = "..."
doi      = "..."
units    = "eV"
orbitals = ["...", ...]   # optional: orbital names in file order (length is checked)

U = [            # Norb x Norb, nested rows (a flat row-major list also works)
  [4.67, 3.52, ...],
  ...
]
J = [ ... ]
```

The bundled sets use the Wannier90 d-orbital order `dz2, dxz, dyz, dx2-y2, dxy` (recorded in
`orbitals`); the orbital order of a file must match that of the Hamiltonian it is used with.
JSON files (`{"U": [...], "J": [...]}`) are still read; a `.toml` of the same name takes precedence.
