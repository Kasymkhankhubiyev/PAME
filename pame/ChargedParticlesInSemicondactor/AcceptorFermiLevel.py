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
from typing import NamedTuple
from pame.ChargedParticlesInSemicondactor.CalculateParticles import *

me_effective = float
mh_effective = float
Kelvin = float
Nc = float
Nd = float
eV = float


class Result(NamedTuple):
    Ef: float
    n: str
    p: str
    Ndneg: str
    Q: float
    ratio: float
    Nv: str
    Nc: str


def find_fermi_level(me: me_effective, mh: mh_effective, Jd: float, t: Kelvin, Efpl: eV, Efneg: eV, Ec: eV, Ev: eV, Na: float):
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

    nc = calc_Nc(me, t)
    nv = calc_Nv(mh, t)

    a, b = Efpl, Efneg

    Ef = (a + b)*0.5

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    naneg = calc_Naneg(Na=Na, Ef=Ef, Ea=Jd - Ev, t=t)
    q = count_Q(n=n, p=p, Na=naneg)

    if np.abs(q/(n + naneg)) < 0.0001:
        # print(f'Na={Na}     nc={nv}')
        return Result(Ef=Ef, n=convert_charges(n), p=convert_charges(p), Ndneg=convert_charges(naneg),
                      Q=q, ratio=(q/(n + naneg)), Nv=convert_charges(nv), Nc=convert_charges(nc))
    else:
        if q < 0:
            return find_fermi_level(me=me, mh=mh, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Na=Na, Efneg=Ef, Efpl=Efpl)
        elif q > 0:
            return find_fermi_level(me=me, mh=mh, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Na=Na, Efneg=Efneg, Efpl=Ef)
