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


class Charge(NamedTuple):
    name: str
    n: int
    p: int
    Nd: int
    Nc: int
    Q: int


class Nparticle(NamedTuple):
    name: str
    body: float
    power: int


def _count_Q(n: int, p: int, Nd: int) -> int:
    """
    :param n: кол-во негативных носителей
    :param p: кол-во положительных носителей
    :param Nd: кол-во ионизированных атомов
    :return: Q - заряд полупроводника
    """
    Q = n + p + Nd
    return Q


def calc_n(nc: Nparticle, Efpl: float, Efneg: float, t: Kelvin) -> Nparticle:
    # k = 1.38 * 10 ** -23  # J/K
    k = 1.38e-16  # эрг/К

    Nc = (nc.body * 10**nc.power)
    expl = np.exp((Efpl - Efneg)/(k * 6.24e11 * t))

    n = Nc * expl

    power = int(np.log10(n))
    body = float(format(n / 10 ** round(np.log10(n)), '.1f'))
    # print(f' n = {n}')
    # print(f'body = {body}')

    return Nparticle(name='n', body=body, power=power)


def calc_p(nd: Nparticle, Efpl: float, Efneg: float, t: Kelvin) -> Nparticle:
    k = 1.38e-16  # эрг/К

    Nd = (nd.body * 10 ** nd.power)
    expl = np.exp((Efpl - Efneg) / (k * 6.24e11 * t))

    p = Nd * expl
    power = int(np.log10(p))
    body = float(format(p / 10 ** round(np.log10(p)), '.1f'))

    return Nparticle(name='p', body=body, power=power)


def _calc_Nc(me: me_effective, t: Kelvin) -> Nparticle:
    k = 1.38 * 10**-23  # J/K
    h = 1.054 * 10**-34  # kg * m /sec^2
    m0 = 9.109 * 10**-31  # kg ~ 0.511MeV

    Nc = 2 * ((2 * np.pi * me * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3

    Nc /= 10**6  # 1/cm^3

    return Nparticle(name='Nc', body=float(format(Nc/ 10 ** round(np.log10(Nc)), '.1f')),
                     power=round(np.log10(Nc)))


def _calc_Nd(mh: mh_effective, t: Kelvin) -> Nparticle:
    k = 1.38 * 10 ** -23  # J/K
    h = 1.054 * 10 ** -34  # kg * m /sec^2
    m0 = 9.109 * 10 ** -31  # kg ~ 0.511MeV

    Nd = 2 * ((2 * np.pi * mh * m0 * k * t)/((2 * np.pi * h)**2)) ** 1.5  # 1/m^3

    Nd /= 10 ** 6  # 1/cm^3

    return Nparticle(name='Nd', body=float(format(Nd / 10 ** round(np.log10(Nc) - 1), '.1f')),
                     power=round(np.log10(Nd)))


def calculate_charges(me: me_effective, mh: mh_effective, t: Kelvin):
    nc = _calc_Nc(me, t)
    # nd = _calc_Nd(mh, t)
    return nc
