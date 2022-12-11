"""
для Si(Eg=1.12eV) p-n переход при температуре T=250K c уровнями легирования Nd=10^16 cm^-3 для n-типа,
Na = 3*10^16 cm^-3 в p-типе определите эффективные плотности состояний Nc, Nv, разницу уровней Ферми в n- и p- типеб
ширину области обеднения W. Эффективная масса электрона в зоне проводимости m_e = 0.36, эффективная масса дырок
в валентной зоне  m_h = 0.81m_0, epsilon = 11.7

Рассчитать вольт-амперную характеристику такого диода, если площадь перехода составляет 1 мм^2. Диффузионные длины
для дырок взять на сайте.

Оценить максимальный ток через такой pn переход, если кристалл расположен в стандартном корпусе ТО-220, который позволяет
рассеивать до 50Вт, после установки на теплоотвод.
"""
import pame.ChargedParticlesInSemicondactor.AcceptorFermiLevel as afl
import pame.ChargedParticlesInSemicondactor.DonorFermiLevel as dfl
from pame.Semiconductors.helper import pn_junction_w_width



def run():
    Si_epsilon = 11.7
    Si_Nd, Si_Na = 1e16, 3e16
    Si_n = dfl.find_fermi_level(me=0.36, mh=0.81, t=250, Jd=0.05, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=Si_Nd)
    print(f'Si_n: Nv={Si_n.Nv}, Nc={Si_n.Nc}')
    Si_p = afl.find_fermi_level(me=0.36, mh=0.81, t=250, Jd=0.05, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Na=Si_Na)
    print(f'Si_p: Nv={Si_p.Nv}, Nc={Si_p.Nc}')
    print(f'Fermi Levels difference: {Si_n.Ef-Si_p.Ef}')

    w, w_p, w_n = pn_junction_w_width(delta_phi=Si_n.Ef-Si_p.Ef, epsilon=Si_epsilon, n0=Si_Nd, p0=Si_Na, explicit=True)

    print(f'w_p = {w_p}, w_n = {w_n}, w = {w}')


