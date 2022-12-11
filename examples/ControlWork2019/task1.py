"""
рассчитайте плотность соединения ZnTe. Параметр решетки ZnTe = 6.1 ангстрема"
"""
from pame.AtomsCalculation.molecule import *

# if __name__ == 'main':
def run():
    lattice = 6.1e-8
    m_a = 65 + 99
    ZnTe = Molecule(a=lattice, m_a=m_a)

    ZnTe_atom_dens = ZnTe.atoms_density()
    ZnTe_atom_num = ZnTe.atoms_number()

    print(ZnTe_atom_dens, ZnTe_atom_num)
