"""
На кремниевой подложке с маркировкой КЭФ-10 сформирован барьер Шоттки с алюминием. Определить ширину ОПЗ и
высоту барьера Шоттки. Работы выхода Никеля составляет 5ю1 эВ, электронное сродство в кремнии 4.05 эВ,
ширина запрещенной зоны 1.12 эВ. Подвижность дырок в кремнии 500 см^2/(B*c)
"""
from pame.FermiLevelPinning.DonorFermiLevel import find_fermi_level
from pame.Semiconductors.helper import w_width, nd_from_mobility, pn_junction_w_width


if __name__ == '__main__':
    Si_me, Si_mh = 0.36, 0.81
    Jd, P_res = 0.05, 10  # Ohm*cm
    Ni_E0, Si_ksi = 5.1, 4.05  # eV
    Si_eg, Si_p_mu = 1.12, 500  # eV, cm^2/V/sec

    Nd = nd_from_mobility(ohm_cm=10, mobility=500)
    print(Nd)

    Si_ef = find_fermi_level(me=Si_me, mh=Si_mh, t=300, Jd=Jd, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=Nd).Ef
    print(Si_ef)

    print(Ni_E0+Si_ef)
    print(Si_ksi + Si_eg)

    delta_E = Ni_E0 + Si_ef - (Si_ksi + Si_eg)  # т.к. у нас Si n-типа
    print(delta_E)
    if delta_E > Si_eg-Si_ef:
        print("Инверсия")
        # в этом случае у нас будто p-n переход - можем найти ширину ОПЗ
        delta_phi = delta_E  # V контактная разность потенциалов
        w = w_width(delta_phi=delta_phi, semicond_epsilon=11.7, carrier=Nd)
        print(f'Ширина изгиба: {w/100} м')
        c = 2.998e10
        volt = 1e8 / c
        print(f'Высота барьера: {delta_phi} eV')
    else:
        print("Обеднение - Уровень Ферми ушел глубоко и носителей заряда стало очень мало")
