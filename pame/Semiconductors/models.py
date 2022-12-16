
class Model(object):
    def __init__(self, lattice, epsilon, electron_affinity, Eg, Jd, Ja, me, mh, me_light, mh_heavy,
                 spin_orbital_splitting, mu_electrons, mu_holes, Dn, Dp, Ln, Lp, Nd, Na):
        """
        :param lattice: lattice constant in cm
        :param epsilon: dielectric constant
        :param electron_affinity: in eV
        :param Eg: energy gap in eV
        :param Jd: ionization energy of donors in eV
        :param Ja: ionization energy of acceptors in eV
        :param me: effective mass of density of states in the conduction band in m0 (free electron mass)
        :param mh: effective mass of density of states in the valence band in m0
        :param me_light: mass of light electron in m0
        :param mh_heavy: mass of heavy hole in m0
        :param spin_orbital_splitting: in eV
        :param mu_electrons: mobility electrons in cm^2 V^-1 sec^-1
        :param mu_holes: mobility holes in cm^2 V^-1 sec^-1
        :param Dn: Diffusion coefficient electrons
        :param Dp: Diffusion coefficient holes
        :param Ln: lifetime length in cm
        :param Lp: lifetime length in cm
        :param Nd: donors concentration in cm^-3
        :param Na: acceptors concentration in cm^-3
        """
        self.lattice, self.epsilon, self.electron_affinity = lattice, epsilon, electron_affinity
        self.spin_orbital_splitting, self.Eg = spin_orbital_splitting, Eg
        self.me, self.mh = me, mh
        self.me_light, self.mh_heavy = me_light, mh_heavy
        self.mu_electrons, self.mu_holes = mu_electrons, mu_holes
        self.Dn, self.Dp, self.Ln, self.Lp = Dn, Dp, Ln, Lp
        self.Jd, self.Ja = Jd, Ja
        self.Nd, self.Na = Nd, Na

    def reset_Lp(self, Lp: float) -> None:
        self.Lp = Lp

    def reset_Ln(self, Ln: float) -> None:
        self.Ln = Ln

    def reset_Nd(self, Nd: float) -> None:
        self.Nd = Nd

    def reset_Na(self, Na: float) -> None:
        self.Na = Na

    def reset_Jd(self, Jd: float) -> None:
        self.Jd = Jd

    def reset_Ja(self, Ja: float) -> None:
        self.Ja = Ja


class Si(Model):
    def __init__(self, Nd=0.0, Na=0.0, Ln=3e-3, Lp=5e-3) -> None:
        """
        :param Ln: with Nd = 1e18 as default
        :param Lp: with Na = 1e18 as default
        :param Na: acceptors concentration in cm^-3
        :param Nd: donors concentration in cm^-3

        Jd: ionization energy of shallow donors of P in eV as default
        Ja: ionization energy of shallow acceptors of B in eV as default
        """
        Model.__init__(self, lattice=5.43e-8, epsilon=11.7, electron_affinity=4.05, Eg=1.12,
                       spin_orbital_splitting=0.044, me=0.35, mh=0.81, me_light=0.19, mh_heavy=0.49,
                       Jd=0.045, Ja=0.045, mu_electrons=1400, mu_holes=450, Dn=36, Dp=12, Ln=Ln, Lp=Lp,
                       Nd=Nd, Na=Na)


class Ge(Model):
    def __init__(self, Ln=0.2, Lp=0.3, Nd=0.0, Na=0.0) -> None:
        """
        :param Ln: pure n-type at 300K as default
        :param Lp: pure p-type at 300k as default
        :param Na: acceptors concentration in cm^-3
        :param Nd: donors concentration in cm^-3

        Jd: ionization energy of shallow donors of P in eV as default
        Ja: ionization energy of shallow acceptors of B in eV as default
        """
        Model.__init__(self, lattice=5.658e-8, epsilon=16.2, electron_affinity=4.0, Eg=0.661,
                       spin_orbital_splitting=0.29, me=0.22, mh=0.34, me_light=0.0815, mh_heavy=0.33,
                       Jd=0.013, Ja=0.011, mu_electrons=3900, mu_holes=1900, Dn=100, Dp=50, Ln=Ln, Lp=Lp,
                       Nd=Nd, Na=Na)


class GaAs(Model):
    def __init__(self, Ln=4e-5, Lp=4e-5, Na=0.0, Nd=0.0) -> None:
        """
        :param Ln: pure n-type material as default
        :param Lp: pure p-type material as default
        :param Na: acceptors concentration in cm^-3
        :param Nd: donors concentration in cm^-3

        Jd: ionization energy of shallow donors of Ge in eV as default
        Ja: ionization energy of shallow acceptors of Zn in eV as default
        """
        Model.__init__(self, lattice=5.65325e-8, epsilon=12.9, electron_affinity=4.07, Eg=1.414,
                       spin_orbital_splitting=0.34, me=0.85, mh=0.53, me_light=0.075, mh_heavy=0.51,
                       Jd=0.006, Ja=0.025, mu_electrons=8500, mu_holes=400, Dn=200, Dp=10, Ln=Ln, Lp=Lp,
                       Nd=Nd, Na=Na)


class GaP(Model):
    def __init__(self, Ln=7e-6, Lp=2e-5, Nd=0.0, Na=0.0) -> None:
        """
        :param Ln: longest lifetime as default
        :param Lp: longest lifetime as default
        :param Na: acceptors concentration in cm^-3
        :param Nd: donors concentration in cm^-3

        Jd: ionization energy of shallow donors of S_p in eV as default
        Ja: ionization energy of shallow acceptors of Zn_Ga in eV as default
        """
        Model.__init__(self, lattice=5.4505e-8, epsilon=11.1, electron_affinity=3.8, Eg=2.26,
                       spin_orbital_splitting=0.08, me=0.79, mh=0.83, me_light=0.22, mh_heavy=0.79,
                       Jd=0.107, Ja=0.0697, mu_electrons=250, mu_holes=150, Dn=6.5, Dp=4, Ln=Ln, Lp=Lp,
                       Nd=Nd, Na=Na)


class InAs(Model):
    def __init__(self, Ln=4.5e-5, Lp=1.5e-5, Nd=0.0, Na=0.0) -> None:
        """
        :param Ln: longest lifetime as default
        :param Lp: longest lifetime as default
        :param Na: acceptors concentration in cm^-3
        :param Nd: donors concentration in cm^-3

        Jd: ionization energy of shallow donors of S, Se, Te, Ge, Si, Sn, Cu in eV as default
        Ja: ionization energy of shallow acceptors of Zn in eV as default
        """
        Model.__init__(self, lattice=6.0583e-8, epsilon=15.15, electron_affinity=4.9, Eg=0.354,
                       spin_orbital_splitting=0.41, me=0.29, mh=0.41, me_light=0.023, mh_heavy=0.026,
                       mu_electrons=4e4, mu_holes=5e2, Jd=1e-3, Ja=0.01, Dn=1e3, Dp=13, Lp=Lp, Ln=Ln,
                       Nd=Nd, Na=Na)


class C(Model):
    """
            :param Ln: longest lifetime as default
            :param Lp: longest lifetime as default
            :param Na: acceptors concentration in cm^-3
            :param Nd: donors concentration in cm^-3

            Jd: Nitrogen most common at the energy level of 1.7eV and 4eV below the bottom of the conduction band
            Ja: Boron is a deep acceptor
            """
    def __init__(self, Ln=1, Lp=1, Nd=0.0, Na=0.0):
        Model.__init__(self, lattice=3.567e-8, epsilon=5.7, electron_affinity=0.0, Eg=5.46,
                       spin_orbital_splitting=6e-3, me=0.57, mh=0.8, me_light=0.36, mh_heavy=2.12,
                       mu_electrons=2200, mu_holes=1800, Dn=57, Dp=46, Ln=Ln, Lp=Lp, Ja=0.37, Jd=1.46,
                       Nd=Nd, Na=Na)


class GaSb(Model):
    def __init__(self, Ln=1e-5, Lp=1e-5, Nd=0.0, Na=0.0) -> None:
        """
        :param Ln: longest lifetime as default
        :param Lp: longest lifetime as default
        :param Na: acceptors concentration in cm^-3
        :param Nd: donors concentration in cm^-3

        Jd: ionization energy of shallow donors of Se in eV as default
        Ja: ionization energy of shallow acceptors of Zn in eV as default
        """
        Model.__init__(self, lattice=6.09593e-8, epsilon=15.7, electron_affinity=4.06, Eg=0.726,
                       spin_orbital_splitting=0.8, me=0.57, mh=0.8, me_light=0.11, mh_heavy=0.4,
                       mu_electrons=3000, mu_holes=1e3, Dn=75, Dp=25, Lp=Lp, Ln=Ln, Jd=0.05, Ja=0.037,
                       Na=Na, Nd=Nd)


class InSb(Model):
    def __init__(self, Ln=1e-5, Lp=1e-5, Nd=0.0, Na=0.0) -> None:
        """
        :param Ln: longest lifetime as default
        :param Lp: longest lifetime as default
        :param Na: acceptors concentration in cm^-3
        :param Nd: donors concentration in cm^-3

        Jd: ionization energy of shallow donors of Se in eV as default
        Ja: ionization energy of shallow acceptors of Zn in eV as default
        """
        Model.__init__(self, lattice=6.479e-8, epsilon=16.8, electron_affinity=4.59, Eg=0.17,
                       spin_orbital_splitting=0.8, me=0.25, mh=0.43, me_light=0.014, mh_heavy=0.43,
                       Jd=7e-4, Ja=0.01, mu_electrons=7.7e4, mu_holes=850, Dn=2e3, Dp=22, Ln=Ln, Lp=Lp,
                       Na=Na, Nd=Nd)
