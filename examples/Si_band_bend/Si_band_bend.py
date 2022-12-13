from pame.BandBend.BandBend import *
from pame.FermiLevelPinning.DonorFermiLevel import find_fermi_level


class Semiconductor:
    def __init__(self, me, mh, Jd, E_g, Nd, Nas, Eas, epsilon, E_f):
        self.me, self.mh, self.Jd, self.E_g = me, mh, Jd, E_g
        self.Nd, self.Nas, self.Eas = Nd, Nas, Eas,
        self.epsilon, self.E_f = epsilon, E_f


if __name__ == '__main__':

    tolerance = 1e-10

    Si = Semiconductor(me=0.36, mh=0.49, Jd=0.045, E_g=1.12, Nd=1e17, Nas=1e8, Eas=0.1, epsilon=11.7, E_f=1.)

    fermi_dichotomy = find_fermi_level(me=Si.me, mh=Si.mh, t=300, Jd=Si.Jd, Efpl=0.57,
                                       Efneg=1.12, Ec=Si.E_g, Nd=Si.Nd, Ev=0)

    print('\nИзгиб зоны с уровнем Ферми полученным методом Дихотомии:\n')
    Si.E_f = fermi_dichotomy[0]

    phi_dichotomy_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                          Ef=Si.E_f, phi0=0, phi1=0.5, method='dichotomy', tolerance=tolerance)

    phi_fixed_point_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                            Ef=Si.E_f, phi0=0.05, method='fixed-point', tolerance=tolerance)

    phi_newtown_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                        Ef=Si.E_f, phi0=0.05, method='newtown', tolerance=tolerance)

    phi_secant_d = calculate_band_bend(epsilon=Si.epsilon, Nd=Si.Nd, t=300, Nas=Si.Nas, Eas=Si.Eas, Eout=0,
                                       Ef=Si.E_f, phi0=0.05, method='secant', tolerance=tolerance, delta=1e-2)