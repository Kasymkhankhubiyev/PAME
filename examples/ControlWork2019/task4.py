"""
Полевом транзисторе Si с p-n переходом (размеры L=10мкм, z=10мкм, a=100мкм), определите напряжение и
ток отсечки V_p, I_p, изобразите семейство ВАХ при V_g=0, V_g = V_p/2 при Nd=10^16 cm^-3
"""
from pame.JFET.calculations import volt_amper_characteristics, Ip


def run():
    Si_mu_e, Si_mu_h = 0.1400, 0.0450  # m^2 V^-1 s^-1
    ip = Ip(z=1e-5, L=1e-5, a=1e-7, mu=Si_mu_e, epsilon=11.7, nd=1e-22)
    volt_amper_characteristics(I_p=ip, path='examples/ControlWork2019/')
