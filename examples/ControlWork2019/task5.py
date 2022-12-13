"""
Рассчитайте разрывы зон и изобразите схематически гетероструктуру ZnSe-GaSb.
ZnSe(a=5.67A, Eg=2.82eV, epsilon=5.9), GaSb(a=6.9A, Eg=0.73eV, epsilon=14.4)
"""
from pame.HeteroStructure.Calculation import *


# if __name__ == 'main':
def run():
    GaSb = SemiCond(name='GaSb', lattice=6.9 * 10**-8, epsilon=14.4, E_g=0.73, spin_orbital_splitting=0.80)
    ZnSe = SemiCond(name='ZnSe', lattice=5.67 * 10**-8, epsilon=5.9, E_g=2.82, spin_orbital_splitting=0.4)

    path = 'examples/ControlWork2019/'
    delta_E_c, delta_E_v = process_heterostructure(wide_band=ZnSe,
                                                   slim_band=GaSb,
                                                   path=path)

