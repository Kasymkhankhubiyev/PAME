from pame.HeteroStructure.Calculation import *
import numpy as np

from fompy.constants import eV, me
from fompy.models import KronigPenneyModel, DiracCombModel
from fompy.units import unit
from math import pi

import matplotlib
import matplotlib.pyplot as plt

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

def run():
    a_AlGaAs = (5.6533 - 0.0078 * 0.4) * 10**-8 # Lattice constant, cm
    epsilon_AlGaAs = 12.9 - 2.84 * 0.4
    E_g_AlGaAs = 1.424 + 1.247 * 0.4  # eV

    a_GaAs = 5.6532e-8  # Lattice constant, cm
    epsilon_GaAs = 12.9
    E_g_GaAs = 1.424  # eV

    w_p_AlGaAs = _plasm_freq(a=a_AlGaAs)
    print(f'Wp for Si:  {w_p_AlGaAs}')
    w_p_GaAs = _plasm_freq(a_GaAs)
    print(f'Wp for Si:  {w_p_GaAs}')

    Eg_AlGaAs_model = _energy(a=a_AlGaAs, epsilon=epsilon_AlGaAs)
    Eg_GaAs_model = _energy(a=a_GaAs, epsilon=epsilon_GaAs)

    print(f'Si Eg = {Eg_AlGaAs_model}')
    print(f'Ge Eg = {Eg_GaAs_model}')

    delta_Eg_model = np.abs(Eg_AlGaAs_model - Eg_GaAs_model) / 2.
    print(f'delta E model = {delta_Eg_model}')

    delta_0_AlGaAs = 0.34 - 0.04 * 0.4
    delta_0_GaAs = 0.34

    print(f'delta_0_Si /3 = {delta_0_AlGaAs / 3}')
    print(f'delta_0_Ge /3 = {delta_0_GaAs / 3}')

    print(f'Si_Eg: {E_g_AlGaAs}')
    print(f'Ge_Eg: {E_g_GaAs}')

    _plot_contact_zone(delta_Eg_model=delta_Eg_model,
                       Eg_outer=E_g_AlGaAs,
                       Eg_inner=E_g_GaAs,
                       delta0_outer=delta_0_AlGaAs,
                       delta0_inner=delta_0_GaAs)

    delta_u = E_g_AlGaAs - E_g_GaAs
    periodic_potential(delta_u)
