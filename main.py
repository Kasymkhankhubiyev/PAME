from task1 import runfile
from ChargedParticlesInSemicondactor.GetFdLevel import *
from ChargedParticlesInSemicondactor.GetFaLevel import *
from SemiconCurrent import CurrentCalculus
import matplotlib.pyplot as plt
from task1 import PeriodicPotential as pp


def calc_Ld()->None:
    # засчет дебаевской длины для Кремния

    epsi0 = 1e-10
    k = 1.38e-023  # J/K
    T = 300.
    e = 1.6e-19  # Кулон
    n0 = 1e17 / 1e6  # 1/м^3

    print((epsi0 * k * T / (e * e * n0)) ** 0.5)


def volt_amper_characteristics() -> None:

    def Id(ud: float, vp: float, ug: float):
        Ip = 0.1  # A ток отсечки
        res = Ip * (3*ud/vp - 2 * ((ud + ug)**1.5)/(vp**1.5) - (ug/vp)**1.5)
        return res

    fig = plt.Figure()
    for i in range(1, 10):
        Id_array, vp_array = [], []
        ud = i
        for j in range(1, 10):
            vp = ud*j/10
            id = Id(ud=ud, vp=vp, ug=1.)
            vp_array.append(vp)
            Id_array.append(id)
        plt.plot(vp_array, Id_array, label=f'ud={ud}')
    plt.legend('best')
    plt.savefig('volt_amper_charact')


def g_m():
    pass


def g_d():
    pass


if __name__ == '__main__':

    # volt_amper_characteristics()
    #
    # # поиск урвня Ферми для доноров
    # calc_donor_fermi_levels()
    #
    # # поиск урвня Ферми для акцепторов
    # calc_acceptor_fermi_levels()

    # рассчет тока в полупроводнике в pn-переходе
    # Nd = 5e17
    # Na = 1e18
    #
    # js = CurrentCalculus.count_Js(t=300., Nd=Nd, Na=Na)
    # print(js)


    # рассчет уровня энергии в периодической яме
    pp.run()



