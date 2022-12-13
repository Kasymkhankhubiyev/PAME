import numpy as np
from pame.exceptions import CantMatchMethod
from fompy.models import KronigPenneyModel
from fompy.util.zeros import find_nth_function_zero

import pame.constants as const

nano = float


class Model(object):

    def __init__(self, m: float, u0: float, a: nano, b: nano) -> None:
        self.m, self.u0 = m * const.m0, u0
        self.a, self.b = a, b

    def _cp(self, energy: float, k: int) -> float:
        """

        :param energy: Energy in Joules
        :param k: a mode number
        :return:
        """
        alpha = np.sqrt(2 * self.m * np.abs(energy) * const.eV_to_J / const.h_bar ** 2)
        betta = np.sqrt(2 * self.m * (energy + self.u0) * const.eV_to_J / const.h_bar ** 2)
        return (np.cos(betta * self.b) * np.cos(alpha * self.a) - np.cos(k * (self.a + self.b))) * \
               (2 * alpha * betta) - (alpha**2 + betta**2) * np.sin(betta * self.b) * np.sin(alpha * self.a)

    def _cp_simply(self, energy: float, k: int) -> float:
        """
        :param energy: Energy in Joules
        :param k: a mode number
        :return:
        """
        alpha = np.sqrt(2 * self.m * np.abs(energy) * const.eV_to_J / const.h_bar ** 2)
        p = self.m * self.u0 * const.eV_to_J * self.b * (self.a + self.b) / const.h_bar ** 2
        summ = self.a + self.b
        return (np.cos(alpha * summ) - np.cos(k * summ)) * ((alpha * summ) + p * np.sin(alpha * summ))

    def _cp_shch(self, energy: float, k: int) -> float:  # with hyperbolistic functions
        """
        :math: $\alpha = \sqrt{\frac{2m|E|}{h^2}}$
        :math: $\beta = \sqrt{\frac{2m(E+U_0}{h^2}}$

        equilibrium condition:

        :math: $\frac{\beta^2 - \alpha^2}{2\alpha\beta}sh(\beta b)sin(\alpha a)+cos(\beta b)sh(\alpha a)-
                        cos(k(a + b)) = 0$
        :param energy: Energy level
        :param k: number of mode
        :return:
        """
        alpha = np.sqrt(2 * self.m * np.abs(energy) * const.eV_to_J / const.h_bar ** 2)
        betta = np.sqrt(2 * self.m * (energy + self.u0) * const.eV_to_J / const.h_bar ** 2)
        return (np.cos(betta * self.b) * np.cosh(alpha * self.a) - np.cos(k * (self.a + self.b))) * \
               (2 * alpha * betta) - (betta**2 - alpha**2) * np.sin(betta * self.b) * np.sinh(alpha * self.a)
        # return (np.cosh(betta * self.b) * np.cos(alpha * self.a) - np.cos(k * (self.a + self.b))) * \
        #        (2 * alpha * betta) - (betta ** 2 - alpha ** 2) * np.sinh(betta * self.b) * np.sin(alpha * self.a)

    def _cp_shch_simpy(self, energy: float, k: int) -> float:
        """
        :param energy:
        :param k:
        :return:
        """
        alpha = np.sqrt(2 * self.m * np.abs(energy) * const.eV_to_J / const.h_bar ** 2)
        summ = self.a + self.b
        p = self.m * self.u0 * const.eV_to_J * self.b * summ / const.h_bar**2
        return (np.cosh(alpha * summ) - np.cos(k * summ)) * (alpha * summ) + p * np.sinh(alpha * summ)

    def _dichotomy(self, f, a, b, k=0, n=300) -> float:
        """
        returns an energy of a free particle on the k_th mode
        :param f: function of linear equilibrium
        :param a: length of an area where u=-u0
        :param b: length of ab area where u=0
        :param k: a mode number
        :param n: iterations
        :return: energy of a free particle on the specific mode
        """
        f_a = f(a, k)
        ab2 = (a + b) / 2
        for i in range(n):
            f_2 = f(ab2, k)
            if f_a * f_2 < 0:
                b = ab2
                ab2 = (a + b) / 2
            else:
                a = ab2
                ab2 = (a + b) / 2
                f_a = f(a, k)
        return ab2

    def calculate_energy_levels(self, method: str, level=0):
        methods = ['classic', 'classic_simple', 'shch', 'shch_simple']
        if method in methods:
            if method == 'classic':
                e_level = self._dichotomy(f=self._cp, a=-0.3, b=0, k=level)
                print(f'Kronig-Penney classic model: {e_level}')
                return e_level
            elif method == 'classic_simple':
                e_level = self._dichotomy(f=self._cp_simply, a=-0.3, b=0, k=level)
                print(f'Kroning-Penney classic simple: {e_level}')
                return e_level
            elif method == 'shch':
                e_level = self._dichotomy(f=self._cp_shch, a=-0.3, b=0, k=level)
                print(f'Kroning-Penney shch model: {e_level}')
                return e_level
            elif method == 'shch_simple':
                e_level = self._dichotomy(f=self._cp_shch_simpy, a=-0.3, b=0, k=level)
                print(f'Kroning-Penney shch simple: {e_level}')
                return e_level
        else:
            raise CantMatchMethod(message=method, methods=methods)