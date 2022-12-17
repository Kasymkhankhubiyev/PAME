"""
Рассчитайте разрывы зон и изобразите схематически гетероструктуру GaP-InP.
GaP(a=5.45A, Eg=2.26eV, epsilon=11.), InP(a=5.86A, Eg=1.34eV, epsilon=12.5)

Рассчитать уровни энергии электронов и тяжелых дырок в сверхрешетке на базе
такой гетероструктуры, если толщина слоев GaP=5нм, а InP=7нм
"""
from pame.HeteroStructure.Calculation import *
from pame.KronigPenney.model import Model
from pame.Semiconductors.helper import laser_lambda


# if __name__ == 'main':
def run():
    GaP_width, InP_width = 5e-9, 7e-9
    InP = SemiCond(name='Inp', lattice=5.86 * 10**-8, epsilon=12.5, E_g=1.34, spin_orbital_splitting=0.11)
    GaP = SemiCond(name='GaP', lattice=5.45 * 10**-8, epsilon=11, E_g=2.26, spin_orbital_splitting=0.08)

    path = 'examples/ControlWork2022/'
    delta_E_c, delta_E_v = process_heterostructure(wide_band=GaP,
                                                   slim_band=InP,
                                                   path=path)

    print(f'delta_Ev = {delta_E_v}')
    print(f'delta_Ec = {delta_E_c}')

    energy_e = Model(m=0.08, u0=delta_E_c, a=GaP_width, b=InP_width).calculate_energy_levels(method='classic_simple')
    heavy_d = Model(m=0.6, u0=delta_E_v, a=GaP_width, b=InP_width).calculate_energy_levels(method='classic_simple')

    print(f'lambda = {laser_lambda(delta_energy=-energy_e + heavy_d)*1e9} нм')
