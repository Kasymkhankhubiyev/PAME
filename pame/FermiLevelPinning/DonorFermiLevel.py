"""
В данной программе производим рассчет
заряженных частицы в полупроводнике.
Точнее считаем Nd-, p, n, Q

Пусть ширина запрещенной зоны Eg = 1.12eV
Jd = Ec - Ed = 50 meV

1) Уравнение электронейтральности:
    n = Nd^+ + p

2) Ef - уровень Ферми, Ec - дно зоны проводимости,
Ev - потолок валентной зоны
    n = Nc * exp(- (Ec - Ef)/kT )
    p = Nv * exp(- (Ef - Ev)/kT )

    Получаем Nd^+ = Nd/ (1 + 4 * exp((Ef - Eg)/kT ))

 Q = Nd^+ + p - n --> заряд в материале
I. Предположим, что Nd_- = 0 -> полупроводник собственный, уровень ферми ушел глубоко:
    Ef+ = (Ec + Ev)/2 + 3/4 * kT * ln(me*/mh*)

II. Уровень Ферми совпадает с границей зоны проводимости при
Nd=10^17 -> статистика не вырожденная:
    Nc(Si, T=300K) ~ 10^19
Т.о. на втором участке заряд гарантировано отрицательный,
    Ef1 = (Ef_+ + Ef_-)/2

"""
from typing import NamedTuple
from pame.FermiLevelPinning.CalculateParticles import *
from pame.exceptions import CantMatchMethod, CatchZeroDelta, CantRunDichotomyMethod, CannotMatchResultFormat

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
    Nv: float
    Nc: float


class Result_xey(NamedTuple):
    Ef: float
    n: str
    p: str
    Ndpl: str
    Q: float
    ratio: float
    Nv: str
    Nc: str


def _dichotomy_method(nc: float, nv: float, t: Kelvin, Jd: eV, Ef_low: eV, Ef_upper: eV, Ec: eV,
                      Ev: eV, Nd: float, tolerance=1e-7):
    """
    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param t: temperature in Kelvin
    :param Jd: ionization energy in eV
    :param Ef_low: lowest Fermi energy level in eV
    :param Ef_upper: upper Fermi energy level in eV
    :param Ec: conduction band energy level in eV
    :param Ev: valence band energy level in eV
    :param Nd: concentration of donors
    :param counter: counts iteration steps
    :param tolerance: acceptable error value
    :return: Fermi level  in eV and iterations steps amount
    """
    a, b = Ef_low, Ef_upper

    Ef = (a + b) / 2.

    n = calc_n(nc=nc, Ef=Ef, Ec=Ec, t=t)
    p = calc_p(nv=nv, Ef=Ef, Ev=Ev, t=t)
    ndpl = calc_Ndplus(Nd=Nd, Ef=Ef, Ed=Ec - Jd, t=t)
    q = count_Q(n=n, p=p, Ndpl=ndpl)

    if np.abs(q / (p + ndpl)) < tolerance:
        return Result(Ef=Ef, n=n, p=p, Ndpl=ndpl, Q=q, ratio=q/(p + ndpl), Nv=nv, Nc=nc)
    else:
        if q > 0:
            return _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Ef_upper=Ef, Ef_low=Ef_low)
        elif q < 0:
            return _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ec=Ec, Ev=Ev, Nd=Nd, Ef_upper=Ef_upper, Ef_low=Ef)


def _fixed_point_method(nc: float, nv: float, nd: float, t: Kelvin, e_d: eV, e_f: eV,
                        e_c: eV, e_v: eV, e_f0: eV, tolerance=1e-7) -> tuple:
    """
    :math: $\phi(x)=f(x)+x$ : $x_{n+1}=x_n - \frac{f(x_n)}{f'(\ksi)}$

    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_d: ionization energy in eV
    :param e_f: current Fermi energy level in eV
    :param e_c: conduction band energy level in eV
    :param e_v: valence band energy level in eV
    :param e_f0: fixed level of Fermi energy in eV
    :param counter: counts iteration steps
    :param tolerance: acceptable error value
    :return: Fermi level  in eV and iterations steps amount
    """
    if np.abs(balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)) < tolerance:
        n = calc_n(nc=nc, Ef=e_f, Ec=e_c, t=t)
        p = calc_p(nv=nv, Ef=e_f, Ev=e_v, t=t)
        ndpl = calc_Ndplus(Nd=nd, Ef=e_f, Ed=e_d, t=t)
        q = count_Q(n=n, p=p, Ndpl=ndpl)
        return Result(Ef=e_f, n=n, p=p, Ndpl=ndpl, Q=q, ratio=q/(p + ndpl), Nv=nv, Nc=nc)
    else:
        diff = diff_balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f0, e_c=e_c, e_v=e_v, e_d=e_d)
        _lambda = 1/diff
        ef_i = -_lambda * balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)
        return _fixed_point_method(nc=nc, nv=nv, nd=nd, t=t, e_d=e_d, e_f=e_f+ef_i, e_c=e_c, e_v=e_v,
                                   e_f0=e_f0, tolerance=tolerance)


def _newtown_method(nc: float, nv: float, nd: float, t: Kelvin, e_d: eV, e_f: eV, e_c: eV, e_v: eV,
                    tolerance=1e-7):
    """
    :math: $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$
    :param nc:  concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_d: ionization energy in eV
    :param e_f: current Fermi energy level in eV
    :param e_c: conduction band energy level in eV
    :param e_v: valence band energy level in eV
    :param counter: counts iteration steps
    :param tolerance: acceptable error value
    :return: Fermi level  in eV and iterations steps amount
    """
    if np.abs(balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)) < tolerance:
        n = calc_n(nc=nc, Ef=e_f, Ec=e_c, t=t)
        p = calc_p(nv=nv, Ef=e_f, Ev=e_v, t=t)
        ndpl = calc_Ndplus(Nd=nd, Ef=e_f, Ed=e_d, t=t)
        q = count_Q(n=n, p=p, Ndpl=ndpl)
        return Result(Ef=e_f, n=n, p=p, Ndpl=ndpl, Q=q, ratio=q/(p + ndpl), Nv=nv, Nc=nc)
    else:
        diff = diff_balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d)
        _lambda = balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_c=e_c, e_v=e_v, e_d=e_d) / diff
        return _newtown_method(nc=nc, nv=nv, nd=nd, t=t, e_d=e_d, e_f=e_f-_lambda, e_c=e_c,
                               e_v=e_v, tolerance=tolerance)


def _secant_method(nc: float, nv: float, nd: float, t: Kelvin, e_d: eV, e_f: eV, e_c: eV, e_v: float,
                   tolerance=1e-7, delta=0.1, e_fi=None):
    """
    :param nc:  concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_d: ionization energy in eV
    :param e_f: fermi level value in eV for {n-1}-th iteration
    :param e_c: conduction band energy level in eV
    :param e_v: valence band energy level in eV
    :param counter: counts iteration steps
    :param tolerance: Fermi level  in eV and iterations steps amount
    :param delta: x_1 - x_0 difference for a first iteration
    :param e_fi: fermi level value in eV for n-th iteration
    :return: Fermi level  in eV and iterations steps amount
    """
    if e_fi is not None:
        if (balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_fi, e_c=e_c, e_v=e_v, e_d=e_d)) < tolerance:
            n = calc_n(nc=nc, Ef=e_f, Ec=e_c, t=t)
            p = calc_p(nv=nv, Ef=e_f, Ev=e_v, t=t)
            ndpl = calc_Ndplus(Nd=nd, Ef=e_f, Ed=e_d, t=t)
            q = count_Q(n=n, p=p, Ndpl=ndpl)
            return Result(Ef=e_f, n=n, p=p, Ndpl=ndpl, Q=q, ratio=q/(p + ndpl), Nv=nv, Nc=nc)
        else:
            factor = delta / (balance_function(nc=nc, nv=nv, nd=nd,
                                               t=t, e_f=e_fi, e_c=e_c,
                                               e_v=e_v, e_d=e_d) - balance_function(nc=nc, nv=nv, nd=nd, t=t,
                                                                                    e_f=e_f, e_c=e_c, e_v=e_v,
                                                                                    e_d=e_d))

            ef_ = e_fi - balance_function(nc=nc, nv=nv, nd=nd, t=t, e_f=e_fi, e_c=e_c, e_v=e_v, e_d=e_d) * factor
            return _secant_method(nc=nc, nv=nv, nd=nd, t=t, e_d=e_d, e_f=e_fi, e_c=e_c, e_v=e_v,
                                  tolerance=tolerance, e_fi=ef_)
    else:
        if delta != 0:
            factor = delta / (balance_function(nc=nc, nv=nv, nd=nd,
                                               t=t, e_f=e_f, e_c=e_c,
                                               e_v=e_v, e_d=e_d) - balance_function(nc=nc, nv=nv, nd=nd,
                                                                                    t=t, e_f=e_f-delta, e_c=e_c,
                                                                                    e_v=e_v, e_d=e_d))
            ef_ = e_f - balance_function(nc=nc, nv=nv, nd=nd,
                                         t=t, e_f=e_f, e_c=e_c,
                                         e_v=e_v, e_d=e_d) * factor
            return _secant_method(nc=nc, nv=nv, nd=nd, t=t, e_f=e_f, e_fi=ef_, e_c=e_c, e_v=e_v, e_d=e_d,
                                  tolerance=tolerance)
        else:
            raise CatchZeroDelta


def fermi_methods() -> list:
    return ['dichotomy', 'newtown', 'fixed-point', 'secant']


def calculate_fermi_level(me: me_effective, mh: mh_effective, t: Kelvin, Jd: eV, Ef0: eV, Ec: eV,
                          Nd: float, result_format='numeric', method='dichotomy', Ev=0., tolerance=1e-7,
                          Ef1=None, delta=1e-3) -> tuple:
    """
    :param result_format:
    :param me: эффективная масса электрона
    :param mh: эффективная масса дырки
    :param t: температура
    :param Jd: энергия ионизации
    :param Ef_low: $E_{f}^{+}$ - нижняя граница
    :param Ef_upper: $E_{f}^{-}$ - верхняя граница
    :param Ec: дно зоны проводимости
    :param Ev: потолок валентной зоны
    :param Nd: концентрация доноров
    :param method: метод поиска уровня ферми
    :return: Fermi level in eV and iterations steps amount
    """
    methods = ['dichotomy', 'newtown', 'fixed-point', 'secant']
    result = None

    nc = calc_Nc(me, t)
    nv = calc_Nv(mh, t)
    try:
        if method == 'dichotomy':
            if Ef1 is not None and Ef1 > Ef0:
                result = _dichotomy_method(nc=nc, nv=nv, t=t, Jd=Jd, Ef_low=Ef0, Ef_upper=Ef1,
                                                         Ec=Ec, Ev=Ev, Nd=Nd, tolerance=tolerance)
            else:
                raise CantRunDichotomyMethod(phi0=Ef0, phi1=Ef1)
        elif method == 'fixed-point':
            result = _fixed_point_method(nc=nc, nv=nv, nd=Nd, t=t, e_d=Ec-Jd, e_f=Ef0,
                                                       e_f0=Ef0, e_c=Ec, e_v=Ev, tolerance=tolerance)
        elif method == 'newtown':
            result = _newtown_method(nc=nc, nv=nv, nd=Nd, t=t, e_d=Ec-Jd,
                                                   e_f=Ef0, e_c=Ec, e_v=Ev, tolerance=tolerance)
        elif method == 'secant':
            result = _secant_method(nc=nc, nv=nv, nd=Nd, t=t, e_d=Ec-Jd, e_f=Ef0, e_c=Ec, e_v=Ev,
                                                  tolerance=tolerance, delta=delta)
        else:
            raise CantMatchMethod(message=method, methods=methods)

        if result_format is 'numeric':
            return result
        elif result_format is 'xey':
            return Result_xey(Ef=result.Ef, n=convert_charges(result.n), p=convert_charges(result.p),
                      Ndpl=convert_charges(result.Ndpl), Q=result.Q, ratio=result.ratio,
                      Nv=convert_charges(nv), Nc=convert_charges(nc))
        else:
            raise CannotMatchResultFormat(available_formats=['numeric', 'xey'])

    except CantMatchMethod as e:
        print(e.args)
    except CantRunDichotomyMethod as e:
        print(e.args)
    except CatchZeroDelta as e:
        print('Delta equals zero')
