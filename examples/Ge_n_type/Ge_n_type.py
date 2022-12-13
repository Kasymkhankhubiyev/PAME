"""
GeSb
"""
from pame.FermiLevelPinning.DonorFermiLevel import *
from pame.Semiconductors.helper import debye_length

Kelvin = float
Nd = 1e17
T = 200.
Eea = 4.0  # eV
Jd = 0.01
Au_E0 = 5.4  # Ev
epsilon0 = 1e-10
Ge_epsilon = 16.2
e = 1.602e-19  # Кулон

if __name__ == '__main__':
    me = 0.22
    mh = 0.34
    Eg = 0.742 - 4.8e-4 * T**2 / (T + 235)  # eV
    print(f'Ширина Запрещенной зоны: {Eg}')
    result = find_fermi_level(me=me, mh=mh, t=T, Efpl=Eg / 2, Jd=Jd, Efneg=Eg, Ec=Eg, Ev=0, Nd=1e17)
    print(f'Уровень Ферми: {result.Ef}')
    print(f'концентрация электронов: {result.n}')
    delta_E = Au_E0 - (Eea + Eg)
    print(delta_E)
    if delta_E - result.Ef > 0:
        print("Инверсия")
        # в этом случае у нас будто p-n переход - можем найти ширину ОПЗ
        delta_phi = 0.5  # V контактная разность потенциалов
        # 1eV = 1,602e-19 J
        w = np.sqrt(2 * delta_phi * epsilon0 * Ge_epsilon / (e * Nd * 1e6))  # перевели в кубометр
        print(f'Ширина изгиба: {w * 1e9} нано метров')
    if delta_E - result.Ef <= 0:
        print("Обеднение - Уровень Ферми ушел глубоко и носителей заряда стало очень мало")

    debay_len = debye_length(Ge_epsilon, Nd, T)
    print(f'Величина дебая: {debay_len*1e6} мкм')  #маикрометры


