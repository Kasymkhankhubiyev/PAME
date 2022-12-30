from pame.HeteroStructure.Calculation import *
from pame.KronigPenney.model import Model


if __name__ == '__main__':
    a_AlGaAs = (5.6533 + 0.0078 * 0.4) * 10**-8 # Lattice constant, cm
    epsilon_AlGaAs = 12.9 - 2.84 * 0.4
    E_g_AlGaAs = 1.424 + 1.247 * 0.4  # eV

    a_GaAs = 5.6532e-8  # Lattice constant, cm
    epsilon_GaAs = 12.9
    E_g_GaAs = 1.424  # eV

    AlGaAs = SemiCond(name='AlGaAs', E_g=E_g_AlGaAs, epsilon=epsilon_AlGaAs,
                      lattice=a_AlGaAs, spin_orbital_splitting=0.34 - 0.04 * 0.4)
    GaAs = SemiCond(name='GaAs', E_g=E_g_GaAs, epsilon=epsilon_GaAs,
                    lattice=a_GaAs, spin_orbital_splitting=0.34)

    delta_u_cond, delta_u_valence = process_heterostructure(wide_band=AlGaAs, slim_band=GaAs)

    print(f'глубина ямы:  {delta_u_cond}')
    print(f'delta in valence: {delta_u_valence}')

    m_free_e_GaAs = 0.063
    m_heavy_hole_GaAs = 0.51

    model = Model(m=m_free_e_GaAs, a=5, b=10, u0=delta_u_cond)
    model.calculate_energy_levels(method='shch')

    model = Model(m=-m_heavy_hole_GaAs, a=5, b=10, u0=-delta_u_valence)
    model.calculate_energy_levels(method='shch')
