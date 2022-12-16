import matplotlib.pyplot as plt
import numpy as np
from typing import NamedTuple
from pame.FermiLevelPinning.CalculateParticles import calc_Nv, calc_Nc, calc_n, calc_p
from pame.constants import k, e

me_effective = float
mh_effective = float
Kelvin = float


class Current(NamedTuple):
    js: float
    jn: float
    jp: float


def count_p_n(ni2: float, nc: float, Ef_n: float, Eg: float, t: float) -> float:
    """
    :math: $p_n0 = \frac{n_i^2}{N_cexp(\frac{E_f_n-E_g}{kT})}$
    :param ni2: electrons concentration in non-degenerate semiconductor
    :param nc: electrons density
    :param Ef_n: fermi energy in n-type in eV
    :param Eg: energy gap in eV
    :param t: temperature in Kelvin
    :return: concentration of protons in n-area
    """
    return ni2 / (nc * np.exp((Ef_n-Eg)/(k * t)))


def count_n_p(ni2: float, nv: float, Ef_p: float, t: float) -> float:
    """
    :math: $n_p0 = \frac{n_i^2}{N_vexp(\frac{-E_f_p}{kT})}$
    :param ni2: electrons concentration in non-degenerate semiconductor
    :param nv: protons density
    :param Ef_p: fermi energy in p-type in eV
    :param t: temperature in eV
    :return: concentration of electrons in p-area
    """
    return ni2 / (nv * np.exp((-Ef_p)/(k * t)))


def count_ni(t: float, me: float, mh: float, Eg: float) -> float:
    Nc = calc_Nc(me=me, t=t)
    Nv = calc_Nv(mh=mh, t=t)
    # k = 1.38e-16  # эрг/К

    return (Nc*Nv * np.exp(-1. * Eg/(2 * k * t)))**0.5


def count_pi(t: float, me: float, mh: float, Eg: float) -> float:
    Nc = calc_Nc(me=me, t=t)
    Nv = calc_Nv(mh=mh, t=t)
    return (Nc * Nv * np.exp(-1. * Eg / (2 * k * t)))**0.5  # 10^3


def _count_Jp(me: float, mh: float, Eg: float, Ef_n: float, t: float, Dp: float, Lp: float) -> float:
    """
    :param me: effective mass of electron
    :param mh: effective mass of proton
    :param Eg: energy gap in eV
    :param Ef_n: fermi energy of n-type in eV
    :param t: temperature in Kelvin
    :param Dp: protons diffusion coefficient
    :param Lp: protons diffusion length
    :return: current density
    """
    pn0 = count_p_n(ni2=count_pi(t=t, me=me, mh=mh, Eg=Eg)**2, nc=calc_Nc(me=me, t=t), Ef_n=Ef_n, Eg=Eg, t=t)
    return (e * Dp * pn0) / Lp


def _count_Jn(me: float, mh: float, Eg: float, Ef_p: float, t: float, Dn: float, Ln: float) -> float:
    """
    :param me: effective mass of electron
    :param mh: effective mass of proton
    :param Eg: energy gap in eV
    :param Ef_p: fermi energy for p-type in eV
    :param t: temperature in Kelvin
    :param Dn: electrons diffusion coefficient
    :param Ln: electrons diffusion length
    :return: current density
    """
    np0 = count_n_p(ni2=count_ni(t=t, me=me, mh=mh, Eg=Eg)**2, nv=calc_Nv(mh=mh, t=t), Ef_p=Ef_p, t=t)

    return (e * Dn * np0) / Ln


def count_Js(me: float, mh: float, t: float, Eg: float,
             Dn: float, Ln: float, Dp: float, Lp: float, Ef_n: float, Ef_p: float) -> Current:

    Jn = _count_Jn(t=t, Dn=Dn, Ln=Ln, me=me, mh=mh, Eg=Eg, Ef_p=Ef_p)
    Jp = _count_Jp(t=t, Dp=Dp, Lp=Lp, me=me, mh=mh, Eg=Eg, Ef_n=Ef_n)
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
    """

    :param js: current density in A/cm^2
    :param s: area in cm^2
    :param t: temperature in Kelvin
    :param path: path to a directory to save a VAC image into
    :return: max current
    """
    k = 1.38e-23  # J/K
    u = np.linspace(-1, 0.1, 1000)
    J = []
    for i in range(len(u)):
        J.append(js * s * (np.exp(u[i]*e/(k*t)) - 1))
    # plt.yscale('log')
    plt.plot(u, J, color='blue')
    plt.savefig(path+'vac')
    return u, J


