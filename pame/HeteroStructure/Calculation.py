import numpy as np
import matplotlib.pyplot as plt
from typing import NamedTuple

"""
Модель кронекера-пелли
сделать периодический потенциал

сверх решетка! 
"""
lattice = float


class SemiCond(NamedTuple):
    name: str
    lattice: float
    epsilon: float
    E_g: float
    spin_orbital_splitting: float


def _plasm_freq(a: lattice) -> float:
    e = 4.803e-10
    m0 = 9.1e-28  # Rest mass of an electron, g
    return np.sqrt(4 * np.pi * (32 / pow(a, 3)) * pow(e, 2) / m0)  # СГС


def _energy(a: lattice, epsilon: float) -> float:
    h = 4.135e-15  # Planck's constant, Ev * c
    erg = 1.6e-12  # 1 eV
    return h * _plasm_freq(a) / np.sqrt(epsilon - 1) * erg / 1e-11


def _plot_contact_zone(delta_Eg_model: float, Eg_outer: float, Eg_inner: float,
                       delta0_outer: float, delta0_inner: float) -> plt:
    x = np.linspace(0, 10, 10)
    plt.plot(x, np.zeros_like(x), color='blue')
    plt.plot(x, np.zeros_like(x) + Eg_outer, color='blue')
    plt.plot(x, np.zeros_like(x) - delta0_outer / 3, 'r--')
    x = np.linspace(10, 20, 10)
    plt.plot(x, np.zeros_like(x) - delta0_outer / 3 + delta_Eg_model, color='green')
    plt.plot(x, np.zeros_like(x) - delta0_outer / 3 + delta_Eg_model + delta0_inner / 3, 'g-.')
    plt.plot(x, np.zeros_like(x) - delta0_outer / 3 + delta_Eg_model + delta0_inner / 3 + Eg_inner, color='green')

    if - delta0_outer / 3 + delta_Eg_model < 0:
        z = np.linspace(- delta0_outer / 3 + delta_Eg_model, 0, 10)
    else:
        z = np.linspace(0, - delta0_outer / 3 + delta_Eg_model, 10)
    plt.plot(np.zeros_like(z) + 10, z, color='orange')

    if - delta0_outer / 3 + delta_Eg_model + delta0_inner / 3 + Eg_inner < Eg_outer:
        z = np.linspace(- delta0_outer / 3 + delta_Eg_model + delta0_inner / 3 + Eg_inner, Eg_outer, 10)
    else:
        z = np.linspace(Eg_outer, - delta0_outer / 3 + delta_Eg_model + delta0_inner / 3 + Eg_inner, 10)
    plt.plot(np.zeros_like(z) + 10, z, color='orange')

    plt.savefig(f'path/model.png')
    return plt


def process_heterostructure(wide_band: SemiCond, slim_band: SemiCond, path=None) -> tuple:
    wide_band_w_p = _plasm_freq(a=wide_band.lattice)
    print(f'Wp for Si:  {wide_band_w_p}')
    slim_band_w_p = _plasm_freq(a=slim_band.lattice)
    print(f'Wp for Si:  {slim_band_w_p}')

    e_g_wide_band_model = _energy(a=wide_band.lattice, epsilon=wide_band.epsilon)
    e_g_slim_band_model = _energy(a=slim_band.lattice, epsilon=wide_band.epsilon)

    print(f'{wide_band.name} Eg = {e_g_wide_band_model}')
    print(f'{slim_band.name} Eg = {e_g_slim_band_model}')

    delta_e_g_model = np.abs(e_g_wide_band_model - e_g_slim_band_model) / 2.
    print(f'delta E model = {delta_e_g_model}')

    print(f'delta_0_{wide_band.name} /3 = {wide_band.spin_orbital_splitting / 3}')
    print(f'delta_0_{slim_band.name} /3 = {slim_band.spin_orbital_splitting / 3}')

    print(f'{wide_band.name}_Eg: {wide_band.E_g}')
    print(f'{wide_band.name}_Eg: {slim_band.E_g}')

    if path is not None:
        plt = _plot_contact_zone(delta_Eg_model=delta_e_g_model,
                                 Eg_outer=wide_band.E_g,
                                 Eg_inner=slim_band.E_g,
                                 delta0_outer=wide_band.spin_orbital_splitting,
                                 delta0_inner=slim_band.spin_orbital_splitting)

    delta_e_conductivity = wide_band.E_g - (0 - wide_band.spin_orbital_splitting / 3 + delta_e_g_model +
                                            slim_band.spin_orbital_splitting / 3 + slim_band.E_g)

    return delta_e_conductivity, delta_e_g_model


def run():
    """
    Валентная зона из 3 частей =
    :return:
    """

    a_Si = 5.431e-8  # Lattice constant, cm
    epsilon_Si = 11.7  # Dielectric constant
    E_g_Si = 1.17  # eV

    a_Ge = 5.6532e-8  # Lattice constant, cm
    epsilon_Ge = 12.9  # Dielectric constant
    E_g_Ge = 0.65  # eV

    w_p_si = _plasm_freq(a=a_Si)
    print(f'Wp for Si:  {w_p_si}')
    w_p_Ge = _plasm_freq(a_Ge)
    print(f'Wp for Si:  {w_p_Ge}')

    Eg_Si_model = _energy(a=a_Si, epsilon=epsilon_Si)
    Eg_Ge_model = _energy(a=a_Ge, epsilon=epsilon_Ge)

    print(f'Si Eg = {Eg_Si_model}')
    print(f'Ge Eg = {Eg_Ge_model}')

    delta_Eg_model = np.abs(Eg_Si_model - Eg_Ge_model) / 2.
    print(f'delta E model = {delta_Eg_model}')

    delta_0_Si = 0.044
    delta_0_Ge = 0.29

    print(f'delta_0_Si /3 = {delta_0_Si / 3}')
    print(f'delta_0_Ge /3 = {delta_0_Ge / 3}')

    print(f'Si_Eg: {E_g_Si}')
    print(f'Ge_Eg: {E_g_Ge}')

    _plot_contact_zone(delta_Eg_model=delta_Eg_model,
                       Eg_outer=E_g_Si,
                       Eg_inner=E_g_Ge,
                       delta0_outer=delta_0_Si,
                       delta0_inner=delta_0_Ge)
