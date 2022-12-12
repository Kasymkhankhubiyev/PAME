import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple
from pame.ChargedParticlesInSemicondactor.CalculateParticles import calc_Nv, calc_Nc, calc_n, calc_p

me_effective = float
mh_effective = float
Kelvin = float

from fompy.units import unit


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


def count_ni(t: float, me: float, mh: float, Eg: float) -> float:

    # Nc = 6.2 * 10**15 * t**1.5
    # Nv = 3.5 * 10**15 * t**1.5
    Nc = calc_Nc(me=me, t=t)
    Nv = calc_Nv(mh=mh, t=t)
    # Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    ni = (Nc*Nv * np.exp(-1. * Eg/(2 * k * 6.24e11 * t)))**0.5

    return ni


def count_pi(t: float, me: float, mh: float, Eg: float) -> float:
    # Nc = 6.2 * 10**15 * (t ** 1.5)
    # Nv = 3.5 * 10**15 * (t ** 1.5)
    Nc = calc_Nc(me=me, t=t)
    Nv = calc_Nv(mh=mh, t=t)
    # Eg = 1.12  # eV
    k = 1.38e-16  # эрг/К

    pi = (Nc * Nv * np.exp(-1. * Eg / (2 * k * 6.24e11 * t)))**0.5  #10^3

    return pi


def _count_Jp(me: float, mh: float, Na: float, Eg: float, Ef_n: float, t: float, Dp: float, Lp: float) -> float:
    e = 1.6e-19  # Кулон
    # Dp = 12  # cm^2/s
    lp = Lp  # 1e-3  # cm
    nn = Na

    pi = count_pi(t=t, me=me, mh=mh, Eg=Eg)

    # pn0 = pi**2 / nn

    pn0 = count_p_n(ni2=pi**2, nc=calc_Nc(me=me, t=t), Ef_n=Ef_n, Eg=Eg, t=t)

    Jp = (e * Dp * pn0) / lp

    return Jp


def _count_Jn(me: float, mh: float, Nd: float, Eg: float, Ef_p: float, t: float, Dn: float, Ln: float) -> float:
    e = 1.6e-19  # Кулон
    # Dn = Dn  # 36  # cm^2/s
    ln = Ln  # 5e-3  # cm
    pp = Nd

    ni = count_ni(t=t, me=me, mh=mh, Eg=Eg)
    # np0 = ni**2 / pp

    np0 = count_n_p(ni2=ni**2, nv=calc_Nv(mh=mh, t=t), Ef_p=Ef_p, t=t)

    Jn = ((e * Dn * np0) / ln)

    return Jn


def count_Js(me: float, mh: float, t: float, Nd: float, Na: float, Eg: float,
             Dn: float, Ln: float, Dp: float, Lp: float, Ef_n: float, Ef_p: float) -> Current:

    Jn = _count_Jn(Nd=Nd, t=t, Dn=Dn, Ln=Ln, me=me, mh=mh, Eg=Eg, Ef_p=Ef_p)
    Jp = _count_Jp(Na=Na, t=t, Dp=Dp, Lp=Lp, me=me, mh=mh, Eg=Eg, Ef_n=Ef_n)
    Js = Jn + Jp
    return Current(js=Js, jn=Jn, jp=Jp)


def calclate_I_amper(js: float, s: float) -> float:
    """
    :math: $I(A)=j_s*S =~ [Кл/сек]$
    :param js: плотность тока
    :param s: площадь в мм
    :return: ток в амперах
    """
    return js * s


def volt_amper_characteristic(js: float, s: float, t: float, path: str):
    e = 1.6e-19  # Кулон
    k = 1.38e-23  # J/K
    u = np.linspace(-1, 1, 1000)
    J = []
    for i in range(len(u)):
        J.append(js * s * (np.exp(u[i]*e/(k*t)) - 1) / 1e4)
    plt.plot(u, J, color='blue')
    plt.savefig(path+'vac')
    return u, J


