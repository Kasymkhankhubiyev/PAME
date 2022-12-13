import matplotlib.pyplot as plt
import numpy as np
# Ip = 0.1  # A ток отсечки


def Ip(z: float, a: float, L: float, mu: float, epsilon: float, nd: float) -> float:
    epsilon0 = 8.8e-14  # F/cm
    e = 1.602e-19
    return z * mu * e**2 * nd**2 * a**3 / (6 * epsilon0 * epsilon * L)


def Vp(a: float, epsilon: float, nd: float) -> float:
    epsilon0 = 8.8e-14  # F/cm
    e = 1.602e-19  # Culon
    return e * nd * a**2 / (2 * epsilon0 * epsilon)


def volt_amper_characteristics(I_p: float, path, ug: np.array, vp: float) -> None:

    def Id(ud: float, vp: float, ug: float):
        res = I_p * (3 * ud / vp - 2 * ((ud + ug) ** 1.5) / (vp ** 1.5) - (ug / vp) ** 1.5)
        return res

    ud = np.linspace(0, vp, 1000)

    for i in range(len(ug)):
        Id_array = []
        for j in range(len(ud)):
            id = Id(ud=ud[j], vp=vp, ug=ug[i])
            Id_array.append(id)
        plt.plot(ud, Id_array, label=f'ud={ud}')
    # plt.legend('best')
    plt.savefig(path + 'volt_amper_charact')


def g_m(Ip: float, ud: float, vp: float, ug: float):
    gm = -Ip * (3 * (ud + ug)**0.5/vp**1.5 - 3 * ug**0.5 / vp**1.5)
    return gm


def g_d(Ip: float, ud: float, vp: float, ug: float):
    gd = Ip * (3/vp - 3 * (ud + ug)**0.5 / vp**1.5)
    return gd
