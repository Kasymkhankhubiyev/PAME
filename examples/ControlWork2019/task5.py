"""
Рассчитайте разрывы зон и изобразите схематически гетероструктуру ZnSe-GaSb.
ZnSe(a=5.67A, Eg=2.82eV, epsilon=5.9), GaSb(a=6.9A, Eg=0.73eV, epsilon=14.4)

Рассчитать уровни энергии электронов и тяжелых дырок в сверхрешетке на базе
такой гетероструктуры, если толщина слоев ZnSe=5нм, а GaSb=7нм
"""
from pame.HeteroStructure.Calculation import *
from pame.KronigPenney.model import Model


# if __name__ == 'main':
def run():
    ZnSe_width, GaSb_width = 5e-9, 7e-9
    GaSb = SemiCond(name='GaSb', lattice=6.9 * 10**-8, epsilon=14.4, E_g=0.73, spin_orbital_splitting=0.80)
    ZnSe = SemiCond(name='ZnSe', lattice=5.67 * 10**-8, epsilon=5.9, E_g=2.82, spin_orbital_splitting=0.80)

    path = 'examples/ControlWork2019/'
    delta_E_c, delta_E_v = process_heterostructure(wide_band=ZnSe,
                                                   slim_band=GaSb,
                                                   path=path)

    energy_e = Model(m=0.041, u0=delta_E_c, a=ZnSe_width, b=GaSb_width).calculate_energy_levels(method='classic_simple')
    heavy_d = Model(m=0.4, u0=delta_E_v, a=ZnSe_width, b=GaSb_width).calculate_energy_levels(method='classic_simple')
