"""
В данной программе производим рассчет
заряженных частицы в полупроводнике.
Точнее считаем Nd-, p, n, Q

Пусть ширина запрещенной зоны Eg = 1.12eV
Jd = Ec - Eg = 50 meV

1) Уравнение электронейтральности:
    n + Nd-= p

2) Ef - уровень Ферми, Ec - дно зоны проводимости,
Ev - потолок валентной зоны
    n = Nc * exp(- (Ec - Ef)/kT )
    p = Nv * exp(- (Ef - Ev)/kT )

    Получаем Nd- = Nd/ (1 + 4 * exp((Ef - Eg)/kT ))

 Q = Nd- + n - p --> заряд в материале
I. Предположим, что Nd- = 0 -> полупроводник собственный:
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
    Ndpl: str
    Q: float
    ratio: float
    Nv: str
    Nc: str


def find_fermi_level(me: me_effective, mh: mh_effective, t: Kelvin, Jd: eV, Efpl: eV, Efneg: eV, Ec: eV, Ev: eV, Nd: float):
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

    #TODO выводить заряды в виде 1e10

    # Jd = 0.05  # eV

    nc = calc_Nc(me, t)
    nv = calc_Nv(mh, t)

    a, b = Efpl, Efneg

    Ef = (a + b) / 2.

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    ndpl = calc_Ndplus(Nd=Nd, Ef=Ef, Ed=Ec - Jd, t=t)
    q = count_Q(n=n, p=p, Nd=ndpl)
    # print(n, p, ndpl)

    if np.abs(q/(p + ndpl)) < 0.0001:
        # print(f'Nd={Nd}     nc={nc}')
        return Result(Ef=Ef, n=convert_charges(n), p=convert_charges(p),
                      Ndpl=convert_charges(ndpl), Q=q, ratio=(q/(p + ndpl)),
                      Nv=convert_charges(nv), Nc=convert_charges(nc))
    else:
        if q > 0:
            return find_fermi_level(me=me, mh=mh, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Ef, Efpl=Efpl)
        elif q < 0:
            return find_fermi_level(me=me, mh=mh, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Efneg=Efneg, Efpl=Ef)
