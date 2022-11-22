import numpy as np
import matplotlib.pyplot as plt

"""
Модель кронекера-пелли
сделать периодический потенциал

сверх решетка! 
"""
lattice = float


def plasm_freq(a: lattice) -> float:
    e = 4.803e-10
    m0 = 9.1e-28  # Rest mass of an electron, g
    return np.sqrt(4 * np.pi * (32 / pow(a, 3)) * pow(e, 2) / m0)  # СГС


def energy(a: lattice, epsilon: float) -> float:
    h = 4.135e-15  # Planck's constant, Ev * c
    erg = 1.6e-12  # 1 eV
    return h * plasm_freq(a) / np.sqrt(epsilon - 1) * erg / 1e-11


def plot_contact_zone(delta_Eg_model: float, Eg_outer: float, Eg_inner: float,
                      delta0_outer: float, delta0_inner: float) -> None:

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
    plt.plot(np.zeros_like(z)+10, z, color='orange')

    if - delta0_outer / 3 + delta_Eg_model + delta0_inner / 3 + Eg_inner < Eg_outer:
        z = np.linspace(- delta0_outer / 3 + delta_Eg_model + delta0_inner / 3 + Eg_inner, Eg_outer, 10)
    else:
        z = np.linspace(Eg_outer, - delta0_outer / 3 + delta_Eg_model + delta0_inner / 3 + Eg_inner, 10)
    plt.plot(np.zeros_like(z) + 10, z, color='orange')

    plt.savefig('Heterojunctions/model.png')


def energy_split_orbital():
    pass

def run():
    """
    Валентная зона из 3 частей =
    :return:
    """


    a_Si = 5.431e-8  # Lattice constant, cm
    epsilon_Si = 11.7  # Dielectric constant
    E_g_Si = 1.17  # eV

    a_Ge = 5.658e-8  # Lattice constant, cm
    epsilon_Ge = 16.2  # Dielectric constant
    E_g_Ge = 0.65  # eV

    w_p_si = plasm_freq(a=a_Si)
    print(f'Wp for Si:  {w_p_si}')
    w_p_Ge = plasm_freq(a_Ge)
    print(f'Wp for Si:  {w_p_Ge}')

    Eg_Si_model = energy(a=a_Si, epsilon=epsilon_Si)
    Eg_Ge_model = energy(a=a_Ge, epsilon=epsilon_Ge)

    print(f'Si Eg = {Eg_Si_model}')
    print(f'Ge Eg = {Eg_Ge_model}')

    delta_Eg_model = np.abs(Eg_Si_model - Eg_Ge_model) / 2.
    print(f'delta E model = {delta_Eg_model}')

    delta_0_Si = 0.044
    delta_0_Ge = 0.29

    print(f'delta_0_Si /3 = {delta_0_Si/3}')
    print(f'delta_0_Ge /3 = {delta_0_Ge/3}')

    print(f'Si_Eg: {E_g_Si}')
    print(f'Ge_Eg: {E_g_Ge}')

    plot_contact_zone(delta_Eg_model=delta_Eg_model,
                      Eg_outer=E_g_Si,
                      Eg_inner=E_g_Ge,
                      delta0_outer=delta_0_Si,
                      delta0_inner=delta_0_Ge)
