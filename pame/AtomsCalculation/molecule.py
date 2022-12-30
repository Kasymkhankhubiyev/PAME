angstrom = float
molar_mass = float


class Molecule:
    """
    :ru:    Рассчитываем плотность атомов в г/см^3 и кол-во атомов в 1 cм^3
            Для рассчета мы используем кубическую гранецентрированную решетку
            типа цинковой обманки. Полагаем, что у нас каждого элемента по 4 атома.

            кол-во атомов рассчитываем так: Восемь вершин, каждый атом в вершине
            является общим для 8 других элементарных решеток, поэтому учитываем эти
            8 атомов с множителем 1/8, так же на каждой гране по 1 атому, но каждая грань
            принадлежит одновременно двум элементарным решеткам, поэтому 6 атомов учитываем
            с множителем 1/2, и внутри решетки есть еще 4 атома, итого:
            8 * 1/8 + 6 * 1/2 + 4 = 8 атомов на одну элементарную решетку, если у нас молекула
            состоит из двух разных молекул, то берем поровну: по 4 каждого элемента.

    :en:    calculates atoms density in g/cm^3 and atoms amount in 1 cm^3
            in calculations we assume four atoms for each element.
    """
    def __init__(self, a: angstrom, m_a: molar_mass) -> None:
        """
        :param a: lattice in cm
        :param m_a: molar mass in g/mole
        """
        self.a, self.m_a, self.v = a, m_a, a**3

    def atoms_density(self) -> float:
        """
        :math: $\ro = \frac{M_A}{N_A}\times\frac{N}{V}$

        :return: atoms density
        """
        n_a = 6.022e23  # 1/mole
        n = 8
        rho = self.m_a * n / (n_a * self.v)  # g/mole / (1/mole * cm^3) = g/cm^3
        return rho

    def atoms_number(self) -> float:
        """
        :math: $N=\ro V\frac{N_A}{M_A}$

        :return: atoms number
        """
        n_a = 6.022e23
        n = self.atoms_density() * self.v * n_a / self.m_a
        return n
