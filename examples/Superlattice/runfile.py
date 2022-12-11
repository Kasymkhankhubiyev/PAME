from pame.HeteroStructure.Calculation import *
import numpy as np

from fompy.constants import eV, me
from fompy.models import KronigPenneyModel, DiracCombModel
from fompy.units import unit
from math import pi

import matplotlib
import matplotlib.pyplot as plt
from pame.HeteroStructure.Calculation import *

matplotlib.rc('axes.formatter', useoffset=False)


def periodic_potential(u):
    m = 0.063 * me
    U0 = -u * eV  # '-' потенцияальная яма
    a = 10 * 1e-9 * unit('nm')
    b = 5 * 1e-9 * unit('nm')

    kp_model = KronigPenneyModel(a, b, U0)
    dc_model = DiracCombModel(a + b, a * U0)

    es = np.linspace(-0.001 * eV, 0.001 * eV, 100000)

    ks = kp_model.get_ks(es, m)  # Array of k * (a+b)
    plt.plot(ks, es / eV, label='Kronig-Penney')

    print(es/eV)

    ks = dc_model.get_ks(es, m)
    plt.plot(ks, es / eV, label='Dirac comb')

    plt.axhline(0, color='k', linestyle='--')
    plt.xlim(0, pi)
    plt.xlabel("k*(a+b)")
    plt.ylabel("Energy")
    plt.legend()
    plt.savefig('Superlattice/potential.png')


if __name__ == '__main__':
    a_AlGaAs = (5.6533 - 0.0078 * 0.4) * 10**-8 # Lattice constant, cm
    epsilon_AlGaAs = 12.9 - 2.84 * 0.4
    E_g_AlGaAs = 1.424 + 1.247 * 0.4  # eV

    a_GaAs = 5.6532e-8  # Lattice constant, cm
    epsilon_GaAs = 12.9
    E_g_GaAs = 1.424  # eV

    AlGaAs = SemiCond(name='AlGaAs', E_g=E_g_AlGaAs, epsilon=epsilon_AlGaAs, lattice=a_AlGaAs)
    GaAs = SemiCond(name='GaAs', E_g=E_g_GaAs, epsilon=epsilon_GaAs, lattice=a_GaAs)

    process_heterostructure(wide_band=AlGaAs, slim_band=GaAs)

    delta_u = E_g_AlGaAs - E_g_GaAs
    periodic_potential(delta_u)
