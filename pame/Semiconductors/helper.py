import numpy as np
import pame.constants as constants


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
    """

    :param delta_phi: band bend in eV
    :param semicond_epsilon: dielectric constant
    :param carrier: carrier concentration in 1/cm^3
    :return: width of SCR in cm

    SCR: Space Charge Region
    """
    return np.sqrt(delta_phi * constants.eV_to_J * 2 * constants.epsilon0 * semicond_epsilon / (constants.e * carrier))


def pn_junction_w_n(delta_phi: float, epsilon: float, n0: float, p0: float):
    """
    :math: $w_n^2=\frac{2\epsilon\epsilon_0\delta\phi p_0}{en_0(n_0+p_0)}$
    :param delta_phi: Ef_p-type - Ef_n-type  [eV]
    :param epsilon: dielectric constant
    :param n0: amount of electrons in a conduction band in n-type
    :param p0: amount of protons on a valence band in a p-type
    :return: width of bend for n-type
    """
    # epsilon0 = 8.8e-14  # F/cm
    # eV = 1.6e-19
    # e = 1.602e-19
    return np.sqrt((2 * constants.epsilon0 * epsilon * p0 * delta_phi * constants.eV_to_J) /
                   (constants.e**2 * n0 * (n0 + p0))) / 100  # m


def pn_junction_w_p(delta_phi: float, epsilon: float, n0: float, p0: float):
    """
        :math: $w_n^2=\frac{2\epsilon\epsilon_0\delta\phi n_0}{ep_0(n_0+p_0)}$
        :param delta_phi: Ef_p-type - Ef_n-type  [eV]
        :param epsilon: dielectric constant
        :param n0: amount of electrons in a conduction band in n-type
        :param p0: amount of protons on a valence band in a p-type
        :return: width of bend for p-type
        """
    # epsilon0 = 8.8e-14  # F/cm
    # e = 1.602e-19
    return np.sqrt((2 * constants.epsilon0 * epsilon * n0 * delta_phi * constants.eV_to_J) /
                   (constants.e**2 * p0 * (n0 + p0))) / 100  # m


def _pn_junction_w_full(delta_phi: float, epsilon: float, n0: float, p0: float) -> float:
    # epsilon0 = 8.8e-14  # F/cm
    # e = 1.602e-19
    return np.sqrt(2 * epsilon * constants.epsilon0 * delta_phi * constants.eV_to_J * (n0 + p0) /
                   (constants.e**2 * n0 * p0)) / 100  # m


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


def hydrogen_model_energy(epsilon: float, m: float) -> float:
    # return e ** 4 * m / (2 * h_bar ** 2) / eps ** 2
    pass


def laser_lambda(delta_energy: float) -> float:
    """
    :math: $\lambda=\frac{hc}{\delta E}
    :param delta_energy: difference between lowest electron energy level and highest proton level in eV
    :return: wave length
    """
    return constants.h_bar * constants.c / (delta_energy * constants.eV_to_J) / 100
