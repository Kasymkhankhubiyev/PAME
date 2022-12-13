import numpy as np
from pame.exceptions import CantMatchMethod
from fompy.models import KronigPenneyModel
from fompy.util.zeros import find_nth_function_zero



class Model(object):

    def __init__(self, m, u0, a, b) -> None:
        self.m_e, self.eV = 9.1e-31, 1.6e-19
        self.h = 6.626e-34  # J/sec
        self.h_bar = self.h / (2 * np.pi)
        self.m, self.u0 = m * self.m_e, u0
        self.a, self.b = a * 1e-9, b * 1e-9

    def _cp(self, energy: float, k: int) -> float:
        """

        :param energy: Energy in Joules
        :param k: a mode number
        :return:
        """
        alpha = np.sqrt(2 * self.m * np.abs(energy) * self.eV / self.h_bar ** 2)
        betta = np.sqrt(2 * self.m * (energy + self.u0) * self.eV / self.h_bar ** 2)
        return (np.cos(betta * self.b) * np.cos(alpha * self.a) - np.cos(k * (self.a + self.b))) * \
               (2 * alpha * betta) - (alpha**2 + betta**2) * np.sin(betta * self.b) * np.sin(alpha * self.a)

    def _cp_simply(self, energy: float, k: int) -> float:
        """
        :param energy: Energy in Joules
        :param k: a mode number
        :return:
        """
        alpha = np.sqrt(2 * self.m * np.abs(energy) * self.eV / self.h_bar ** 2)
        p = self.m * self.u0 * self.eV * self.b * (self.a + self.b) / self.h_bar ** 2
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
        alpha = np.sqrt(2 * self.m * np.abs(energy) * self.eV / self.h_bar ** 2)
        betta = np.sqrt(2 * self.m * (energy + self.u0) * self.eV / self.h_bar ** 2)
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
        alpha = np.sqrt(2 * self.m * np.abs(energy) * self.eV / self.h_bar ** 2)
        summ = self.a + self.b
        p = self.m * self.u0 * self.eV * self.b * summ / self.h_bar**2
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

    def find_energy_level(self):
        @np.vectorize
        def get_energy(k):
            try:
                return model.get_energy(k, self.m, (e_start, e_end), accuracy)
            except ValueError:
                return None

        m = 0.49 * self.m_e
        a = 5 * 1e-9
        b = 10 * 1e-9
        u0 = 0.58 * self.eV

        model = KronigPenneyModel(a, b, u0)

        e_0 = model.u0*(1 - 1/100)
        accuracy = abs(model.u0)/10000000

        e_start, e_end = self.find_band_range(model, m, e_0)
        print("E_min, eV = ", e_start / self.eV)
        print("E_max, eV = ", e_end / self.eV)

        # kas = np.linspace(-pi, pi, 1000)
        # es = get_energy(kas / model.period)
        # kas, es = filter_gaps(kas, es)

    def find_band_range(self, model, m, e0):
        xtol_coarse = abs(model.u0)/1e5
        xtol_fine = abs(model.u0)/1e7

        def f(energy):
            return model.equation_left_part(energy, m) - 1
        e_start = find_nth_function_zero(f, e0, xtol_coarse, xtol_fine, num=-1)

        def f(energy):
            return model.equation_left_part(energy, m) + 1
        e_end = find_nth_function_zero(f, e_start, xtol_coarse, xtol_fine, num=0)

        return e_start, e_end

    def filter_gaps(self, ks, es):
        _es, _ks = [], []
        for e, k in zip(es, ks):
            if e is not None:
                _es.append(e)
                _ks.append(k)
        return np.array(_ks), np.array(_es)
