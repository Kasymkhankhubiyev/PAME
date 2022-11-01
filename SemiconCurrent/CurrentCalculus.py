import numpy as np
from typing import NamedTuple


class Current(NamedTuple):
    js: float
    jn: float
    jp: float


def count_ni(t: float) -> float:

    Nc = 6.2 * 10**15 * t**1.5
    Nv = 3.5 * 10**15 * t**1.5
    Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    ni = (Nc*Nv)**0.5 * np.exp(-1. * Eg/(k * 6.24e11 * t))

    return ni


def count_pi(t: float) -> float:
    Nc = 6.2 * 10**15 * (t ** 1.5)
    Nv = 3.5 * 10**15 * (t ** 1.5)
    Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    pi = (Nc * Nv) ** 1. * np.exp(-1. * Eg / (k * 6.24e11 * t))  #10^3

    return pi


def count_Jp(Na: float, t: float) -> float:
    e = 1.6e-19  # Кулон
    Dp = 12  # cm^2/s
    lp = 1e-3  # cm
    nn = Na

    pi = count_pi(t=t)

    pn0 = pi**2 / nn

    pn0 = 200

    Jp = (e * Dp * pn0) / lp

    return Jp


def count_Jn(Nd: float, t: float) -> float:
    e = 1.6e-19  # Кулон
    Dn = 36  # cm^2/s
    ln = 5e-3  # cm
    pp = Nd

    ni = count_pi(t=t)
    np0 = ni**2 / pp

    np0 = 100

    Jn = (e * Dn * np0) / ln

    return Jn


def count_Js(t: float, Nd: float, Na: float) -> Current:
    Jn, Jp = count_Jn(Nd=Nd, t=t), count_Jp(Na=Na, t=t)

    Js = Jn + Jp

    return Current(js=Js, jn=Jn, jp=Jp)
