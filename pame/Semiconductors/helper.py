import numpy as np


def debye_length(epsilon: float, n: float, t: float) -> float:
    """
    :param epsilon: dielectric constant
    :param n: carriers concentration
    :param t: temperature in Kelvin
    :return: Debey length in cm
    """
    k = 1.381e-16  # arg/K
    # e = 1  # eV
    e = 4.803e-10

    return np.sqrt(epsilon * k * 6.24e11 * t / (4 * np.pi * e ** 2 * n))


def w_width(delta_phi: float, semicond_epsilon: float, carrier: float) -> float:
    epsilon0 = 1e-10
    # e = 1  # eV  1.602e-19  # Кулон
    e = 4.803e-10
    return np.sqrt(semicond_epsilon * delta_phi / (2 * np.pi * e * carrier))
    # return np.sqrt(delta_phi * 2 * epsilon0 * semicond_epsilon / (e * carrier))  # перевели в кубометр


def pn_junction_w_n(delta_phi: float, epsilon: float, n0: float, p0: float):
    """
    :math: $w_n^2=\frac{2\epsilon\epsilon_0\delta\phi p_0}{en_0(n_0+p_0)}$
    :param delta_phi: Ef_p-type - Ef_n-type  [eV]
    :param epsilon: dielectric constant
    :param n0: amount of electrons in a conduction band in n-type
    :param p0: amount of protons on a valence band in a p-type
    :return: width of bend for n-type
    """
    epsilon0 = 8.8e-14  # F/cm
    e = 1.602e-19
    return np.sqrt((2 * epsilon0 * epsilon * p0 * delta_phi) / (e * n0 * (n0 + p0))) / 100  # m


def pn_junction_w_p(delta_phi: float, epsilon: float, n0: float, p0: float):
    """
        :math: $w_n^2=\frac{2\epsilon\epsilon_0\delta\phi n_0}{ep_0(n_0+p_0)}$
        :param delta_phi: Ef_p-type - Ef_n-type  [eV]
        :param epsilon: dielectric constant
        :param n0: amount of electrons in a conduction band in n-type
        :param p0: amount of protons on a valence band in a p-type
        :return: width of bend for p-type
        """
    epsilon0 = 8.8e-14  # F/cm
    e = 1.602e-19
    return np.sqrt((2 * epsilon0 * epsilon * n0 * delta_phi) / (e * p0 * (n0 + p0))) / 100  # m


def _pn_junction_w_full(delta_phi: float, epsilon: float, n0: float, p0: float) -> float:
    epsilon0 = 8.8e-14  # F/cm
    e = 1.602e-19
    return np.sqrt(2 * epsilon * epsilon0 * delta_phi * (n0 + p0) / (e * n0 * p0)) / 100  # m


def pn_junction_w_width(delta_phi: float, epsilon: float, n0: float, p0: float, explicit=False):
    if explicit:
        w_n = pn_junction_w_n(delta_phi=delta_phi, epsilon=epsilon, n0=n0, p0=p0)
        w_p = pn_junction_w_p(delta_phi=delta_phi, epsilon=epsilon, n0=n0, p0=p0)
        w = w_p + w_n
        return w, w_p, w_n
    else:
        return _pn_junction_w_full(delta_phi=delta_phi, epsilon=epsilon, n0=n0, p0=p0)


def lambda_electron(me: float, t: float):
    return 132 * np.sqrt(1/me) / t**0.5

def periodic_potential_solver():
    pass
