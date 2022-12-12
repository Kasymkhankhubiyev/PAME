import numpy as np
from typing import NamedTuple
from pame.ChargedParticlesInSemicondactor.CalculateParticles import calc_Nv, calc_Nc

me_effective = float
mh_effective = float
Kelvin = float


class Current(NamedTuple):
    js: float
    jn: float
    jp: float


def count_p_n(ni2: float, nc: float, Ef_n: float, Eg: float, t: float) -> float:
    k = 1.38e-16  # эрг/К
    return ni2 / (nc * np.exp(Ef_n-Eg)/(k * 6.24e11 * t))


def count_n_p(ni2: float, nv: float, Ef_p: float, t: float) -> float:
    k = 1.38e-16  # эрг/К
    return ni2 / (nv * np.exp(-Ef_p)/(k * 6.24e11 * t))


def count_ni(t: float, me: float, mh: float) -> float:

    # Nc = 6.2 * 10**15 * t**1.5
    # Nv = 3.5 * 10**15 * t**1.5
    Nc = calc_Nc(me=me, t=t)
    Nv = calc_Nv(mh=mh, t=t)
    Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    ni = (Nc*Nv)**0.5 * np.exp(-1. * Eg/(2 * k * 6.24e11 * t))

    return ni


def count_pi(t: float, me: float, mh: float) -> float:
    # Nc = 6.2 * 10**15 * (t ** 1.5)
    # Nv = 3.5 * 10**15 * (t ** 1.5)
    Nc = calc_Nc(me=me, t=t)
    Nv = calc_Nv(mh=mh, t=t)
    Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    pi = (Nc * Nv) ** .5 * np.exp(-1. * Eg / (2 * k * 6.24e11 * t))  #10^3

    return pi


def _count_Jp(me: float, mh: float, Na: float, t: float, Dp: float, Lp: float) -> float:
    e = 1.6e-19  # Кулон
    # Dp = 12  # cm^2/s
    lp = Lp  # 1e-3  # cm
    nn = Na

    pi = count_pi(t=t, me=me, mh=mh)

    pn0 = pi**2 / nn

    # pn0 = 200

    Jp = (e * Dp * pn0) / lp

    return Jp


def _count_Jn(me: float, mh: float, Nd: float, t: float, Dn: float, Ln: float) -> float:
    e = 1.6e-19  # Кулон
    # Dn = Dn  # 36  # cm^2/s
    ln = Ln  # 5e-3  # cm
    pp = Nd

    ni = count_pi(t=t, me=me, mh=mh)
    np0 = ni**2 / pp

    # np0 = 100

    Jn = (e * Dn * np0) / ln

    return Jn


def count_Js(me: float, mh: float, t: float, Nd: float, Na: float, Dn: float, Ln: float, Dp: float, Lp: float) -> Current:
    Jn = _count_Jn(Nd=Nd, t=t, Dn=Dn, Ln=Ln, me=me, mh=mh)
    Jp = _count_Jp(Na=Na, t=t, Dp=Dp, Lp=Lp, me=me, mh=mh)
    Js = Jn + Jp
    return Current(js=Js, jn=Jn, jp=Jp)
