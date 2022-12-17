import matplotlib.pyplot as plt
import numpy as np

from pame.constants import epsilon0, e


def Ip(z: float, a: float, L: float, mu: float, epsilon: float, nd: float) -> float:
    """
    :math: $\frac{z\mu e^2N_d^2a^3}{6*\epsilon\epsilon_0L}$

    :param z: depth in cm
    :param a: width in cm
    :param L: channel length in cm
    :param mu: electrons mobility cm^2/(V*sec)
    :param epsilon: dielectric constant
    :param nd: atoms concentration
    :return: ток насыщения
    """
    return z * mu * e**2 * nd**2 * a**3 / (6 * epsilon0 * epsilon * L)


def Vp(a: float, epsilon: float, nd: float) -> float:
    """
    :math: $\frac{eN_d^2a^2}{2\epsilon\epsilon_0}$
    :param a:
    :param epsilon:
    :param nd:
    :return:
    """
    return e * nd * a**2 / (2 * epsilon0 * epsilon)


def volt_amper_characteristics(I_p: float, path, ug: np.array, vp: float) -> None:

    def Id(ud: float, vp: float, ug: float):
        return I_p * (3 * ud / vp - 2 * ((ud + ug) ** 1.5) / (vp ** 1.5) - (ug / vp) ** 1.5)

    ud = np.linspace(0, vp*1.2, 100)

    for i in range(len(ug)):
        Id_array = []
        for j in range(len(ud)):
            # проверяем если $U_d+U_g > U_p$, ток не растет!
            if ud[j] + ug[i] <= vp:
                id = Id(ud=ud[j], vp=vp, ug=ug[i])
                Id_array.append(id)
            else:
                Id_array.append(Id_array[j-1])
        print(f'start current Vg={ug[i]}  Id = {Id_array[0]}')
        print(f'stop current Vg ={ug[i]}, Id = {Id_array[90]}')
        plt.plot(ud, Id_array, label=f'Ug={ug[i]}')
    plt.xlabel('Ug[V]')
    plt.ylabel('Id[A]')
    plt.legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')
    plt.savefig(path + 'volt_amper_charact')


def g_m(Ip: float, ud: float, vp: float, ug: float):
    gm = -Ip * (3 * (ud + ug)**0.5/vp**1.5 - 3 * ug**0.5 / vp**1.5)
    return gm


def g_d(Ip: float, ud: float, vp: float, ug: float):
    gd = Ip * (3/vp - 3 * (ud + ug)**0.5 / vp**1.5)
    return gd
