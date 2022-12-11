import numpy as np


def debye_length(epsilon: float, n: float, t: float) -> float:
    """
    :param epsilon: dielectric constant
    :param n: carriers concentration
    :param t: temperature in Kelvin
    :return: Debey length in cm
    """
    k = 1.381e-16  # arg/K
    e = 1  # eV

    return np.sqrt(epsilon * k * 6.24e11 * t / (4 * np.pi * e ** 2 * n))


def periodic_potential_solver():
    pass
