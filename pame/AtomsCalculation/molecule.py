angstrom = float
molar_mass = float


class Molecule:
    def __init__(self, a: angstrom, m_a: molar_mass) -> None:
        self.a, self.m_a, self.v = a, m_a, a**3

    def atoms_density(self) -> float:
        n_a = 6.022e23
        n = 8
        rho = self.m_a * n / (n_a * self.v)
        return rho

    def atoms_number(self) -> float:
        n_a = 6.022e23
        n = self.atoms_density() * self.v * n_a / self.m_a
        return n
