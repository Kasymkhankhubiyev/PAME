from pame.HeteroStructure.Calculation import *


# creating components of a heterostructure
# Lattice constant in cm
# Energy gap in eV
Si = SemiCond(name='Si', lattice=5.431e-8,
              E_g=1.12, epsilon=11.7,
              spin_orbital_splitting=0.044)

Ge = SemiCond(name='Si', lattice=5.6532e-8,
              E_g=0.65, epsilon=12.9,
              spin_orbital_splitting=0.29)

path = 'path_to_directory_to_save_a_fig_into'
delta_e_c, delta_e_v = process_heterostructure(wide_band=Si, slim_band=Ge, path=path)
