from pame.HeteroStructure.Calculation import *
import numpy as np

from fompy.constants import eV, me
from fompy.models import KronigPenneyModel, DiracCombModel
from fompy.units import unit
from math import pi

import matplotlib
import matplotlib.pyplot as plt
from pame.HeteroStructure.Calculation import *
from pame.KronigPenney.model import Model

matplotlib.rc('axes.formatter', useoffset=False)


def periodic_potential(u):
    m = 0.85 * me
    U0 = -u * eV  # '-' потенцияальная яма
    a = 10 * unit('nm')
    b = 5 * unit('nm')

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
    plt.savefig('examples/Superlattice/potential.png')


# if __name__ == '__main__':
def run():
    a_AlGaAs = (5.6533 + 0.0078 * 0.4) * 10**-8 # Lattice constant, cm
    epsilon_AlGaAs = 12.9 - 2.84 * 0.4
    E_g_AlGaAs = 1.424 + 1.247 * 0.4  # eV

    a_GaAs = 5.6532e-8  # Lattice constant, cm
    epsilon_GaAs = 12.9
    E_g_GaAs = 1.424  # eV

    AlGaAs = SemiCond(name='AlGaAs', E_g=E_g_AlGaAs, epsilon=epsilon_AlGaAs,
                      lattice=a_AlGaAs, spin_orbital_splitting=0.34 - 0.04 * 0.4)
    GaAs = SemiCond(name='GaAs', E_g=E_g_GaAs, epsilon=epsilon_GaAs,
                    lattice=a_GaAs, spin_orbital_splitting=0.34)

    delta_u_cond, delta_u_valence = process_heterostructure(wide_band=AlGaAs, slim_band=GaAs)

    print(f'глубина ямы:  {delta_u_cond}')
    print(f'delta in valence: {delta_u_valence}')

    m_free_e_GaAs = 0.063
    m_heavy_hole_GaAs = 0.51

    model = Model(m=m_free_e_GaAs, a=5, b=10, u0=delta_u_cond)
    model.calculate_energy_levels(method='shch')

    model = Model(m=-m_heavy_hole_GaAs, a=5, b=10, u0=-delta_u_valence)
    model.calculate_energy_levels(method='shch')

    # model.find_energy_level()
