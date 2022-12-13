"""
На кремниевой подложке с маркировкой КДБ-5 сформирован барьер Шоттки с алюминием. Определить ширину ОПЗ и
высоту барьера Шоттки. Работы выхода алюминия составляет 4.2 эВ, электронное сродство в кремнии 4.05 эВ,
ширина запрещенной зоны 1.12 эВ. Подвижность дырок в кремнии 500 см^2/(B*c)
"""
from pame.FermiLevelPinning.AcceptorFermiLevel import find_fermi_level
from pame.Semiconductors.helper import w_width

from fompy.models import *
from fompy.materials import Si
from fompy.units import unit, from_unit
from fompy.constants import *


# if __name__ == 'main':
def run():
    Si_me, Si_mh = 0.36, 0.81
    Jd, B_res = 0.045, 5  # Ohm*cm
    Al_E0, Si_ksi = 4.2, 4.05  # eV
    Si_eg, Si_p_mu = 1.12, 500  # eV, cm^2/V/sec

    qulon = 2997924579.99957

    Nd = 1 / (5 * 500 * e) * qulon
    print(Nd)

    # sem = DopedSemiconductor.from_materials(Si, from_unit(450, 'cm / V s'), 'B', from_unit(1, 'Ohm cm'))
    # print(f'Na: {sem.Na:e}, Es: {sem.Ea / eV}, Ef: {sem.fermi_level() / eV}')

    Si_ef = find_fermi_level(me=Si_me, mh=Si_mh, t=300, Jd=Jd, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Na=Nd).Ef
    print(Si_ef)

    delta_E = (Si_ksi + Si_eg) - Al_E0  # т.к. у нас Si n-типа
    if delta_E - Si_ef > 0:
        print("Инверсия")
        # в этом случае у нас будто p-n переход - можем найти ширину ОПЗ
        delta_phi = delta_E  # V контактная разность потенциалов
        # 1eV = 1,602e-19 J
        w = w_width(delta_phi=delta_phi, semicond_epsilon=11.7, carrier=Nd)
        # w = np.sqrt(2 * delta_phi * epsilon0 * Ge_epsilon / (e * Nd * 1e6))  # перевели в кубометр
        print(f'Ширина изгиба: {w}')
        c = 2.998e10
        volt = 1e8 / c
        print(f'Высота барьера: {delta_phi} eV и {delta_phi/volt} V')
    if delta_E - Si_ef <= 0:
        print("Обеднение - Уровень Ферми ушел глубоко и носителей заряда стало очень мало")
