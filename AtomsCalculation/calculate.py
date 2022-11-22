class Molecule:
    def __init__(self, a, m_a) -> None:
        self.a, self.m_a, self.v = a, m_a, a**3


def atoms_density(molecule: Molecule):
    n_a = 6.022e23
    n = 8
    rho = molecule.m_a * n / (n_a * molecule.v)
    return rho


def atoms_number(molecule: Molecule):
    n_a = 6.022e23
    n = atoms_density(molecule=molecule) * molecule.v * n_a / molecule.m_a
    return n
