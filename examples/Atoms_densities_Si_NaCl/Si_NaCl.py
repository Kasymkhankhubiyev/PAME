from pame.AtomsCalculation.molecule import Molecule

# create a Silicium Molecule
# Lattice in angstrom and mass is a molar mass
Si = Molecule(a=5.43, m_a=28.085)

Si_atom_dens = Si.atoms_density()
Si_atom_num = Si.atoms_number()


# create a NaCl Molecule
NaCl = Molecule(a=5.64, m_a=58.5)

NaCl_atom_dens = NaCl.atoms_density()
NaCl_atom_num = NaCl.atoms_number()
