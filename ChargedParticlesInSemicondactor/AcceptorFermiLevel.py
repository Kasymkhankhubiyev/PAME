"""
В данной программе производим рассчет
заряженных частицы в полупроводнике.
Точнее считаем Nd+, p, n, Q

Пусть ширина запрещенной зоны Eg = 1.12eV
Jd = Ec - Eg = 50 meV

1) Уравнение электронейтральности:
    n = Nd+ + p

2) Ef - уровень Ферми, Ec - дно зоны проводимости,
Ev - потолок валентной зоны
    n = Nc * exp(- (Ec - Ef)/kT )
    p = Nv * exp(- (Ef - Ev)/kT )

    Получаем Nd+ = Nd/ (1 + 1/2 * exp((Ef - Ea)/kT ))

 Q = n - p - Nd+ - заряд в материале
I. Предположим, что Nd+ = 0 -> полупроводник собственный:
    Ef+ = (Ec + Ev)/2 + 3/4 * kT * ln(me*/mh*)

II. Уровень Ферми совпадает с границей зоны проводимости при
Nd=10^17 -> статистика не вырожденная:
    Nc(Si, T=300K) ~ 10^19
Т.о. на втором участке заряд гарантировано отрицательный,
    Ef1 = (Ef+ - Ef-)/2

"""

import numpy as np
from typing import NamedTuple

me_effective = float
mh_effective = float
Kelvin = float
Nc = float
Nd = float
eV = float


class Result(NamedTuple):
    Ef: float
    n: float
    p: float
    Ndneg: float
    Q: float
    ratio: float


def _count_Q(n: float, p: float, Na: float) -> float:
    """
    :param n: кол-во негативных носителей
    :param p: кол-во положительных носителей
    :param Nd: кол-во ионизированных атомов
    :return: Q - заряд полупроводника
    """

    Q = p - n - Na
    return Q


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


def calc_Naneg(Na: float, Ef: eV, Ea: eV, t: Kelvin):
    k = 1.38e-16  # эрг/К

    naneg = Na / (1. + 4. * np.exp((Ef - Ea)/(k * 6.24e11 * t)))
    return naneg


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


def calculate_charges(me: me_effective, mh: mh_effective, t: Kelvin, Efpl: eV, Efneg: eV, Ec: eV, Ev: eV, Nd: float):
    """
    Расчет делаем методом дихотомии

    :param me: эффективная масса плотности состояний электронов
    :param mh: эффективная масса плотности состояний дырок
    :param t: температура
    :param Efpl: уровень Ферми совпадает с дном зоны проводимости,
    все электроны в зоне проводимости, а дырки в валентной
    :param Efng: уровень Ферми близок к потолку валентной зоны
    :return:
    """
    Jd = 0.05

    nc = _calc_Nc(me, t)
    nv = _calc_Nv(mh, t)

    a, b = Efpl, Efneg

    Ef = (a + b)*0.5

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    naneg = calc_Naneg(Na=Nd, Ef=Ef, Ea=Jd-Ev, t=t)
    q = _count_Q(n=n, p=p, Na=naneg)


    if np.abs(q/(n + naneg)) < 0.0001:
        print(f'Na={Nd}     nc={nv}')
        return Result(Ef=Ef, n=n, p=p, Ndneg=naneg, Q=q, ratio=(q/(n + naneg)))
    else:
        if q < 0:
            return calculate_charges(me=me, mh=mh, t=t, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Ef, Efpl=Efpl)
        elif q > 0:
            return calculate_charges(me=me, mh=mh, t=t, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Efneg, Efpl=Ef)
