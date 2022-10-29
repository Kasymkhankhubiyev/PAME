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

    Получаем Nd+ = Nd/ (1 + 1/2 * exp((Ef - Eg)/kT ))

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
    Ndpl: float
    Q: float
    ratio: float


def _count_Q(n: float, p: float, Nd: float) -> float:
    """
    :param n: кол-во негативных носителей
    :param p: кол-во положительных носителей
    :param Nd: кол-во ионизированных атомов
    :return: Q - заряд полупроводника
    """
    Q = n - p - Nd
    return Q


def calc_n(nc: float, Ef: float, Ec: float, t: Kelvin) -> float:  # Nparticle:
    """
    n = Nc * exp(- (Ec - Ef)/kT)
    """
    k = 1.38e-16  # эрг/К

    Nc = nc  # (nc.body * 10**nc.power)
    expl = np.exp((Ef - Ec)/(k * 6.24e11 * t))

    n = Nc * expl

    power = int(np.log10(n))
    body = float(format(n / 10 ** round(np.log10(n)), '.1f'))

    return n    # Nparticle(name='n', body=body, power=power)


def calc_p(nv: float, Ef: float, Ev: float, t: Kelvin) -> float:  # Nparticle:
    """
    p = Nv * exp((Ev - Ef)/kT)
    """
    k = 1.38e-16  # эрг/К

    Nv = nv  #(nv.body * 10 ** nv.power)
    expl = np.exp((Ev - Ef) / (k * 6.24e11 * t))

    p = Nv * expl

    power = int(np.log10(p))
    body = float(format(p / 10 ** round(np.log10(p)), '.1f'))

    return p  # Nparticle(name='p', body=body, power=power)

def calc_Ndplus(Nd: float, Ef: eV, Eg: eV, t: Kelvin):
    k = 1.38e-16  # эрг/К

    ndpl = Nd / (1. + 0.5 * np.exp((Ef - Eg)/ (k * 6.24e11 * t)))

    return ndpl


def _calc_Nc(me: me_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38 * 10**-23  # J/K
    h = 1.054 * 10**-34  # kg * m /sec^2
    m0 = 9.109 * 10**-31  # kg ~ 0.511MeV

    Nc = 2 * ((2 * np.pi * me * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3

    Nc /= 10**6  # 1/cm^3

    return Nc
    # Nparticle(name='Nc', body=float(format(Nc/ 10 ** round(np.log10(Nc)), '.1f')), power=round(np.log10(Nc)))


def _calc_Nv(mh: mh_effective, t: Kelvin) -> float:  # Nparticle:
    k = 1.38 * 10 ** -23  # J/K
    h = 1.054 * 10 ** -34  # kg * m /sec^2
    m0 = 9.109 * 10 ** -31  # kg ~ 0.511MeV

    Nv = 2 * ((2 * np.pi * mh * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3

    Nv /= 10 ** 6  # 1/cm^3

    return Nv
    # Nparticle(name='Nv', body=float(format(Nv / 10 ** round(np.log10(Nv) - 1), '.1f')), power=round(np.log10(Nv)))


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

    nc = _calc_Nc(me, t)
    nv = _calc_Nv(mh, t)

    a, b = Efpl, Efneg

    Ef = (a + b) / 2.

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    ndpl = calc_Ndplus(Nd=Nd, Ef=Ef, Eg=Ec, t=t)
    q = _count_Q(n=n, p=p, Nd=ndpl)


    if np.abs(q/(p + ndpl)) < 0.001:
        return Result(Ef=Ef, n=n, p=p, Ndpl=ndpl, Q=q, ratio=(q/(p + ndpl)))
    else:
        if q > 0:
            return calculate_charges(me=me, mh=mh, t=t, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Ef, Efpl=Efpl)
        elif q < 0:
            return calculate_charges(me=me, mh=mh, t=t, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Efneg, Efpl=Ef)
