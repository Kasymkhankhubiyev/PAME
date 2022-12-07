from pame.HeteroStructure.Calculation import *

# creating components of a heterostructure
# Lattice constant in cm
# Energy gap in eV
AlGaAs = SemiCond(name='AlGaAs',
                  lattice=(5.6533 - 0.0078 * 0.4) * 10**-8,
                  epsilon=12.9 - 2.84 * 0.4,
                  spin_orbital_splitting=0.34 - 0.04 * 0.4,
                  E_g=1.424 + 1.247 * 0.4)

GaAs = SemiCond(name='GaAs',
                lattice=5.6532e-8,
                epsilon=12.9,
                E_g=1.424,
                spin_orbital_splitting=0.34)
path = 'path_to_directory_to_save_a_fig_into'
delta_E_c, delta_E_v = process_heterostructure(wide_band=AlGaAs, slim_band=GaAs, path=path)
