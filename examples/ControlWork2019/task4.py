"""
Полевом транзисторе Si с p-n переходом (размеры L=10мкм, z=10мкм, a=100мкм), определите напряжение и
ток отсечки V_p, I_p, изобразите семейство ВАХ при V_g=0, V_g = V_p/2 при Nd=10^16 cm^-3

10мкм = 1000мксм = 1000 / 10**6 = 1е-3
100мкм = 1е4мксм = 1е4/1е6 = 1е-2
"""
import numpy as np
from pame.JFET.calculations import volt_amper_characteristics, Ip, Vp


def run():
    Si_mu_e, Si_mu_h = 0.1400, 0.0450  # m^2 V^-1 s^-1
    ip = Ip(z=1e-3, L=1e-3, a=1e-2, mu=Si_mu_e, epsilon=11.7, nd=1e16)
    vp = Vp(a=1e-2, epsilon=11.7, nd=1e16)
    print(f'I_p = {ip}, V_p = {vp}')
    volt_amper_characteristics(I_p=ip, path='examples/ControlWork2019/', vp=vp, ug=np.array([0, vp/2]))