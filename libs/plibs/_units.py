#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Physical constants used by the transport / quantum-oscillation drivers, with a switch to a
dimensionless test mode in which every constant is 1.
"""
from typing import NamedTuple

import scipy.constants as scconst

# Public API of this submodule; keeps `from ._x import *` in __init__.py from
# re-exporting numpy, scipy, ctypes and the sibling packages imported above.
__all__ = [
    'Units', 'unit_constants'
]


class Units(NamedTuple):
    """Constant set handed to the drivers (all 1 in the dimensionless mode)."""
    hbar: float      #: reduced Planck constant [eV s]
    ihbar: float     #: 1e-10/hbar: v = ihbar * dE/dk with E in eV, k in AA^-1 -> m/s
    kb: float        #: Boltzmann constant [eV/K]
    eC: float        #: elementary charge [C]
    tau_unit: float  #: unit of the relaxation time: 1 fs = 1e-15 s
    emass: float     #: electron rest mass [kg]
    si: bool         #: True for SI/eV units, False for the dimensionless test mode

    @property
    def tau_label(self) -> str:
        """Unit string printed after a relaxation time ('fs', or '' when dimensionless)."""
        return 'fs' if self.si else ''


def unit_constants(sw_unit: bool = True) -> Units:
    """
    @brief Physical constants in SI/eV units, or all set to 1 for dimensionless code checks.
    @param sw_unit: True for physical units, False for the dimensionless test mode
    @return units: Units(hbar, ihbar, kb, eC, tau_unit, emass, si)
    """
    if not sw_unit:
        return Units(1., 1., 1., 1., 1., 1., False)
    hbar = scconst.physical_constants['reduced Planck constant in eV s'][0]
    return Units(hbar=hbar,
                 ihbar=1.0e-10/hbar,
                 kb=scconst.physical_constants['Boltzmann constant in eV/K'][0],
                 eC=scconst.e,
                 tau_unit=1.e-15,
                 emass=scconst.m_e,
                 si=True)
