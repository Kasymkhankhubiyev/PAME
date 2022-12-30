"""
для Si(Eg=1.12eV) p-n переход при температуре T=200K c уровнями легирования Nd=10^17 cm^-3 для n-типа,
Na = 2*10^16 cm^-3 в p-типе определите эффективные плотности состояний Nc, Nv, разницу уровней Ферми в n- и p- типеб
ширину области обеднения W. Эффективная масса электрона в зоне проводимости m_e = 0.36, эффективная масса дырок
в валентной зоне  m_h = 0.81m_0, epsilon = 11.7

Рассчитать вольт-амперную характеристику такого диода, если площадь перехода составляет 1 мм^2. Диффузионные длины
для дырок взять на сайте.

Оценить максимальный ток через такой pn переход, если кристалл расположен в стандартном корпусе ТО-220, который позволяет
рассеивать до 50Вт, после установки на теплоотвод.

из приближения водородноподобной модели находим энергии Еа и Ed
"""
import pame.FermiLevelPinning.AcceptorFermiLevel as afl
import pame.FermiLevelPinning.DonorFermiLevel as dfl
from pame.Semiconductors.helper import pn_junction_w_width, hydrogen_model_energy, w_width
from pame.SemiconCurrent.CurrentCalculus import count_Js, volt_amper_characteristic
from pame.FermiLevelPinning.CalculateParticles import calc_n, calc_p, calc_Nc, calc_Nv


if __name__ == '__main__':
    Si_epsilon = 11.7
    Si_Nd, Si_Na = 1e17, 2e16
    Si_Dp, Si_Lp, Si_Dn, Si_Ln = 12, 8e-3, 36, 5e-2

    Ea = hydrogen_model_energy(epsilon=Si_epsilon, m=0.81)
    Jd = hydrogen_model_energy(epsilon=Si_epsilon, m=0.36)
    Ed = 1.12 - hydrogen_model_energy(epsilon=Si_epsilon, m=0.36)
    print(f'Ea = {Ea}')
    print(f'Jd = {Jd}')
    print(f'Ed = {Ed}')

    Si_n = dfl.find_fermi_level(me=0.36, mh=0.81, t=200, Jd=Jd, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=Si_Nd)
    print(f'Si_n: Nv={Si_n.Nv}, Nc={Si_n.Nc}')
    Si_p = afl.find_fermi_level(me=0.36, mh=0.81, t=200, Jd=Ea, Efpl=0.0, Efneg=0.58, Ec=1.12, Ev=0, Na=Si_Na)
    print(f'Si_p: Nv={Si_p.Nv}, Nc={Si_p.Nc}')
    print(f'Si_n fermi: {Si_n.Ef}, Si_p fermi: {Si_p.Ef}')
    print(f'Fermi Levels difference: {Si_n.Ef-Si_p.Ef}')

    print(f'Si_n: n0 = {Si_n.n}, Si_p: p0 = {Si_p.p}')

    w, w_p, w_n = pn_junction_w_width(delta_phi=Si_n.Ef - Si_p.Ef, epsilon=Si_epsilon,
                                      n0=calc_n(nc=calc_Nc(me=0.36, t=200), Ef=Si_n.Ef, Ec=1.12, t=200),
                                      p0=calc_p(nv=calc_Nv(mh=0.81, t=200), Ef=Si_p.Ef, Ev=0, t=200),
                                      explicit=True)
    print(f'w_p = {w_p}, w_n = {w_n}, w = {w}')

    # print(f'W = {w_width(delta_phi=Si_n.Ef - Si_p.Ef, semicond_epsilon=Si_epsilon, )}')

    j_current = count_Js(me=0.36, mh=0.81, t=200, Eg=1.12, Dn=Si_Dn, Ln=Si_Ln, Dp=Si_Dp, Lp=Si_Lp, Ef_n=Si_n.Ef,
                         Ef_p=Si_p.Ef)
    print(f'Js = {j_current.js}, j_p = {j_current.jp}, j_n = {j_current.jn}')

    I = volt_amper_characteristic(js=j_current.jn, s=1e-2, path='examples/ControlWork2022/', t=200)
    # print(f'I = {I} Ампер')
