import numpy as np

me_effective = float
mh_effective = float
Kelvin = float
Nc = float
Nd = float
eV = float


def count_Q(n: float, p: float, Ndpl=None, Naneg=None) -> float:
    """
    :param n: кол-во негативных носителей
    :param p: кол-во положительных носителей
    :param Ndpl: кол-во ионизированных атомов
    :return: Q - заряд полупроводника
    """

    if Ndpl is not None:
        return n - p - Ndpl
    if Naneg is not None:
        return n + Naneg - p


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


def calc_Ndplus(Nd: float, Ef: eV, Ed: eV, t: Kelvin):
    k = 1.38e-16  # эрг/К

    ndpl = Nd / (1. + 0.5 * np.exp((Ef - Ed) / (k * 6.24e11 * t)))
    return ndpl


def calc_Naneg(Na: float, Ef: eV, Ea: eV, t: Kelvin):
    k = 1.38e-16  # эрг/К
    naneg = Na / (1. + .25 * np.exp((Ea - Ef)/(k * 6.24e11 * t)))
    return np.abs(naneg)


def calc_Nc(me: me_effective, t: Kelvin) -> float:  # Nparticle:
    """
    :math: $N_c=2(\frac{2\pi m_ckT}{h^2}^{3/2})$
    :param me: effective mass of electron
    :param t: temperature in Kelvin
    :return: electrons concentration 1/cm^3
    """
    k = 1.38e-023  # J/K
    h = 1.054e-034  # J*sec ~~ kg * m / sec
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nc = 2 * ((2 * np.pi * me * m0 * k * t) / (h ** 2)) ** 1.5  # 1/m^3
    Nc /= 10**6  # 1/cm^3
    return Nc


def calc_Nv(mh: mh_effective, t: Kelvin) -> float:  # Nparticle:
    """
    :math: $N_v=2(\frac{2\pi m_ckT}{h^2}^{3/2})$
    :param me: effective mass of electron
    :param t: temperature in Kelvin
    :return: protons concentration 1/cm^3
    """
    k = 1.38e-023  # J/K
    h = 1.054e-034  # kg * m /sec^2
    m0 = 9.109e-031  # kg ~ 0.511MeV

    Nv = 2 * ((2 * np.pi * mh * m0 * k * t) / (h ** 2)) ** 1.5  # 1/m^3
    Nv /= 10 ** 6  # 1/cm^3
    return Nv


def convert_charges(charge: float) -> str:
    factor = round(np.log10(charge))
    body = charge / 10 ** factor
    return '{:.2f}'.format(body)+f'e{factor}'


def count_nc_nv(m_eff: float, t: float) -> float:
    return 2.51e19 * m_eff**1.5 * (t/300)**1.5


def balance_function(nc: float, nv: float, nd: float, t: Kelvin, e_f: eV, e_c: eV, e_v: eV, e_d: eV) -> float:
    """
    :math: $Q=n - (N_d^+ + p)$
    :math: $Q = N_c \factor exp(\frac{E_f - E_c}{kT}) - N_v \factor exp(\frac{E_v-E_f}{kT}) - N_d \factor
            \frac{1}{1+0.5exp(\frac{E_f - E_d}{kT})}$
    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_f: Fermi energy level in eV
    :param e_c: Conduction band energy level in eV
    :param e_v: Valence band energy level in eV
    :param e_d: Ionization energy in eV
    :return: difference in positively and negatively charged particles number
    """
    k = 1.38e-16  # эрг/К

    n = nc * np.exp((e_f - e_c) / (k * 6.24e11 * t))
    p = nv * np.exp((e_v - e_f) / (k * 6.24e11 * t))
    nd_plus = nd / (1. + 0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t)))
    Q = n - nd_plus - p
    return Q/(p + nd_plus)


def diff_balance_function(nc: float, nv: float, nd: float, t: Kelvin, e_f: eV, e_c: eV, e_v: eV, e_d: eV) -> float:
    """
    $\frac{Q}{p + N_d^+}$
    :param nc: concentration of electrons
    :param nv: concentration of holes
    :param nd: concentration of donors
    :param t: temperature in Kelvin
    :param e_f: Fermi energy level in eV
    :param e_c: Conduction band energy level in eV
    :param e_v: Valence band energy level in eV
    :param e_d: Ionized Donors energy level in eV
    :return: differential of difference in positively and
             negatively charged particles number for a specific fermi level
    """
    k = 1.38e-16  # эрг/К
    n = nc * np.exp((e_f - e_c) / (k * 6.24e11 * t))
    p = nv * np.exp((e_v - e_f) / (k * 6.24e11 * t))
    nd_plus = nd / (1. + 0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t)))
    dn = nc * np.exp((e_f - e_c) / (k * 6.24e11 * t)) / (k * 6.24e11 * t)
    dp = -nv * np.exp((e_v - e_f) / (k * 6.24e11 * t)) / (k * 6.24e11 * t)
    dnd_plus = - nd / (1. + 0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t))) ** 2 * \
               0.5 * np.exp((e_f - e_d) / (k * 6.24e11 * t)) / (k * 6.24e11 * t)

    return (dn*(p + nd_plus) - (dp + dnd_plus)*n)/(p + nd_plus)**2
