from task1 import runfile
from ChargedParticlesInSemicondactor.DonorFermiLevel import *
# from ChargedParticlesInSemicondactor.AcceptorFermiLevel import *
from SemiconCurrent import CurrentCalculus


def calc_Ld()->None:
    # засчет дебаевской длины для Кремния

    epsi0 = 1e-10
    k = 1.38e-023  # J/K
    T = 300.
    e = 1.6e-19  # Кулон
    n0 = 1e17 / 1e6  # 1/м^3

    print((epsi0 * k * T / (e * e * n0)) ** 0.5)


def calc_fermi_levels() -> None:
    # поиск урвня Ферми для доноров
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=1e16)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=2 * 1e16)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=4 * 1e16)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=8 * 1e16)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=1e17)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=2 * 1e17)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=4 * 1e17)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=8 * 1e17)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=1e18)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=2 * 1e18)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=4 * 1e18)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=8 * 1e18)
    print(result)
    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=1e19)
    print(result)



if __name__ == '__main__':

    # поиск урвня Ферми для доноров
    calc_fermi_levels()

    # поиск урвня Ферми для акцепторов
    # result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=10e17)
    # print(result)

    # рассчет тока в полупроводнике в pn-переходе
    # Nd = 5e17
    # Na = 1e18
    #
    # js = CurrentCalculus.count_Js(t=300., Nd=Nd, Na=Na)
    # print(js)

    # рассчет уровня энергии в периодической яме
    # runfile.run()



