"""
Рассчитайте плотность антимонида индия InSb. Параметр решетки антимонида индия 6.479A
"""
from pame.AtomsCalculation.molecule import *

def run():
    lattice = 6.479e-8  # cm
    molar = 236.578 / 2  # g/mole
    InSb = Molecule(a=lattice, m_a=molar)

    InSb_atom_dens = InSb.atoms_density()

    print(f'density = {InSb_atom_dens}')