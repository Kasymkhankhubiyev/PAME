import numpy as np

me_effective = float
mh_effective = float
Kelvin = float
Nc = float
Nd = float
eV = float


def count_Q(n: float, p: float, Nd=None, Na=None) -> float:
    """
    :param n: кол-во негативных носителей
    :param p: кол-во положительных носителей
    :param Nd: кол-во ионизированных атомов
    :return: Q - заряд полупроводника
    """

    if Nd is not None:
        return n - p - Nd
    if Na is not None:
        return p - n - Na


def calc_n(nc: float, Ef: float, Ec: float, t: Kelvin) -> float:  # Nparticle:
    """
    n = Nc * exp(- (Ec - Ef)/kT)
    """
    k = 1.38e-16  # эрг/К

    expl = np.exp((Ef - Ec)/(k * 6.24e11 * t))
    n = nc * expl
    return n


def calc_p(nv: float, Ef: float, Ev: float, t: Kelvin) -> float:  # Nparticle:
    """
    p = Nv * exp((Ev - Ef)/kT)
    """
    k = 1.38e-16  # эрг/К

    expl = np.exp((Ev - Ef) / (k * 6.24e11 * t))
    p = nv * expl
    return p


def calc_Ndplus(Nd: float, Ef: eV, Ed: eV, t: Kelvin):
    k = 1.38e-16  # эрг/К

    ndpl = Nd / (1. + 0.5 * np.exp((Ef - Ed) / (k * 6.24e11 * t)))
    return ndpl


def calc_Naneg(Na: float, Ef: eV, Ea: eV, t: Kelvin):
    k = 1.38e-16  # эрг/К

    naneg = Na / (1. + 4. * np.exp((Ef - Ea)/(k * 6.24e11 * t)))
    return naneg


def calc_Nc(me: me_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38e-023  # J/K
    h = 1.054e-034  # kg * m /sec^2
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nc = 2 * ((2 * np.pi * me * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3
    Nc /= 10**6  # 1/cm^3
    return Nc


def calc_Nv(mh: mh_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38e-023  # J/K
    h = 1.054e-034  # kg * m /sec^2
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nv = 2 * ((2 * np.pi * mh * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3
    Nv /= 10 ** 6  # 1/cm^3
    return Nv


def convert_charges(charge: float) -> str:
    factor = round(np.log10(charge))
    body = charge / 10 ** factor
    # print(body, factor, charge)
    # print('{:.2f}'.format(body)+f'e{factor}')
    return '{:.2f}'.format(body)+f'e{factor}'
