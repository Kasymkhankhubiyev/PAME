import numpy as np
from typing import NamedTuple

me_effective = float
mh_effective = float
Kelvin = float


class Current(NamedTuple):
    js: float
    jn: float
    jp: float


def _calc_Nc(me: me_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38e-023  # J/K
    h = 1.054e-034  # kg * m /sec^2
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nc = 2 * ((2 * np.pi * me * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3
    Nc /= 10**6  # 1/cm^3
    return Nc


def _calc_Nv(mh: mh_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38e-023  # J/K
    h = 1.054e-034  # kg * m /sec^2
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nv = 2 * ((2 * np.pi * mh * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3
    Nv /= 10 ** 6  # 1/cm^3
    return Nv


def _count_ni(t: float) -> float:

    me_si = 0.36
    mh_si = 0.81

    # Nc = 6.2 * 10**15 * t**1.5
    # Nv = 3.5 * 10**15 * t**1.5
    Nc = _calc_Nc(me=me_si, t=t)
    Nv = _calc_Nv(mh=mh_si, t=t)
    Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    ni = (Nc*Nv)**0.5 * np.exp(-1. * Eg/(2 * k * 6.24e11 * t))

    return ni


def _count_pi(t: float) -> float:
    me_si = 0.36
    mh_si = 0.81
    # Nc = 6.2 * 10**15 * (t ** 1.5)
    # Nv = 3.5 * 10**15 * (t ** 1.5)
    Nc = _calc_Nc(me=me_si, t=t)
    Nv = _calc_Nv(mh=mh_si, t=t)
    Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    pi = (Nc * Nv) ** .5 * np.exp(-1. * Eg / (2 * k * 6.24e11 * t))  #10^3

    return pi


def _count_Jp(Na: float, t: float) -> float:
    e = 1.6e-19  # Кулон
    Dp = 12  # cm^2/s
    lp = 1e-3  # cm
    nn = Na

    pi = _count_pi(t=t)

    pn0 = pi**2 / nn

    # pn0 = 200

    Jp = (e * Dp * pn0) / lp

    return Jp


def _count_Jn(Nd: float, t: float) -> float:
    e = 1.6e-19  # Кулон
    Dn = 36  # cm^2/s
    ln = 5e-3  # cm
    pp = Nd

    ni = _count_pi(t=t)
    np0 = ni**2 / pp

    # np0 = 100

    Jn = (e * Dn * np0) / ln

    return Jn


def count_Js(t: float, Nd: float, Na: float) -> Current:
    Jn, Jp = _count_Jn(Nd=Nd, t=t), _count_Jp(Na=Na, t=t)

    Js = Jn + Jp

    return Current(js=Js, jn=Jn, jp=Jp)
