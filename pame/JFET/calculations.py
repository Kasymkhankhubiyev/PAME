import matplotlib.pyplot as plt
import numpy as np

from pame.constants import epsilon0, e


def Ip(z: float, a: float, L: float, mu: float, epsilon: float, nd: float) -> float:
    """
    :ru: Ток насыщения
    :en: Saturation current

    :math: $\frac{z\mu e^2N_d^2a^3}{6*\epsilon\epsilon_0L}$

    :param z: transistor depth in cm
    :param a: transistor width in cm
    :param L: transistor channel length in cm
    :param mu: electrons mobility cm^2/(V*sec)
    :param epsilon: dielectric constant
    :param nd: atoms concentration
    :return: ток насыщения
    """
    return z * mu * e**2 * nd**2 * a**3 / (6 * epsilon0 * epsilon * L)


def Vp(a: float, epsilon: float, nd: float) -> float:
    """
    :ru: Напряжение отсечки
    :en: Saturation voltage

    :math: $\frac{eN_d^2a^2}{2\epsilon\epsilon_0}$

    :param a: width in cm
    :param epsilon:
    :param nd:
    :return: Saturation voltage in Volts
    """
    return e * nd * a**2 / (2 * epsilon0 * epsilon)


def Id(I_p: float, ud: float, vp: float, ug: float) -> float:
    """
    :ru: Ток сток-исток
    :en: Drain-Source current

    :math: $I_d = I_p\times(\frac{3U_d}{V_p} - 2(\frac{U_d+U_g}{})^{3/2} - (\frac{U_g}{V_p})^{3/2})$

    :param ud: drain voltage in Volts
    :param vp: saturation voltage in Volts
    :param ug: gate voltage in Volts
    :return: drain-source current in Amper
    """
    return I_p * (3 * ud / vp - 2 * ((ud + ug) ** 1.5) / (vp ** 1.5) - (ug / vp) ** 1.5)


def volt_amper_characteristics(I_p: float, ug: np.array, vp: float, path=None, ud_steps=100) -> None:
    """
    :ru:    Вольт-Амперная характеристика полевого транзистора.

            Ток будет расти до тек пор, пока суммарное напряжение $U_d + U_g$ не достигнет
            напряжения насыщения $V_p$. После этого с ростом напряжения ток остает константным, т.е.
            попадаем в область насыщения.

    :en:    Volt-Amper characteristics of a given JFET

            Current grows until total voltage magnitude $U_d + U_g$ reaches $V_p$,
            after that with voltage growth the current remains constant,
            representing saturation area

    :param I_p: saturation current in Amper
    :param path: path to a directory to save a VAC into
    :param ug: a numpy array of possible values for a gate voltage in Volts
    :param vp: saturation voltage
    :param: ud_steps: number of steps to divide U_d

    :if the path is not None returns a picture of a VAC
    """
    ud = np.linspace(0, vp*1.2, ud_steps)

    for i in range(len(ug)):
        Id_array = []
        for j in range(len(ud)):
            # проверяем если $U_d+U_g > U_p$, ток не растет!
            if ud[j] + ug[i] <= vp:
                id = Id(I_p=I_p, ud=ud[j], vp=vp, ug=ug[i])
                Id_array.append(id)
            else:
                Id_array.append(Id_array[j-1])
        plt.plot(ud, Id_array, label=f'Ug={ug[i]}')
    plt.xlabel('Ug[V]')
    plt.ylabel('Id[A]')
    plt.legend(fontsize=7,
                  ncol=1,
                  facecolor='oldlace',
                  edgecolor='r')

    if path is not None:
        plt.savefig(path + 'volt_amper_characteristic.png')


def g_m(Ip: float, ud: float, vp: float, ug: float) -> float:
    """
    :ru: крутизна канала
    :en: channel steepness

    :math: $g_m = \frac{\partial I_d}{\partial U_g}$
    :math: $g_m = -3I_p(\frac{\sqrt{U_d+U_g}}{V_p^{3/2}}+\frac{1}{2}\frac{\sqrt{U_g}}{V_p^{3/2}})$

    :param Ip: saturation current in Amper
    :param ud: drain voltage in Volts
    :param vp: saturation voltage in Volts
    :param ug: gate voltage in Volts
    :return: channel steepness
    """
    gm = -Ip * (3 * (ud + ug)**0.5/vp**1.5 + 3/2 * ug**0.5 / vp**1.5)
    return gm


def g_d(Ip: float, ud: float, vp: float, ug: float) -> float:
    """
    :ru: проводимость канала
    :en: transfer conductance

    :math: $g_m = \frac{\partial I_d}{\partial U_d}$
    :math: $g_m = 3I_p(\frac{1}{V_p}-\frac{\sqrt{U_d+U_g}}{V_p^{3/2}})$

    :param Ip: saturation current in Amper
    :param ud: drain voltage in Volts
    :param vp: saturation voltage in Volts
    :param ug: gate voltage in Volts
    :return: transfer conductance in Ohm^-1
    """
    gd = Ip * (3/vp - 3 * (ud + ug)**0.5 / vp**1.5)
    return gd
