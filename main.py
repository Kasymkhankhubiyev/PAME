from task1 import runfile
from ChargedParticlesInSemicondactor.chargedparticles import *

if __name__ == '__main__':

    # runfile.run()
    nc = calculate_charges(me=0.36, t=300, mh=0.81)
    # nd = calc_n()
    print(nc)

    pl = (0.57 + 1.12)/2.

    Efleft = 0.57
    Efright = 1.12

    step = (Efright - Efleft) / 4.0

    # TODO сделать runfile, в который просто все передадим, а тот все соберет
    """
    Нужно просто разбить линейно область множества энергий Ферми
    передаем в расчет точку и считаем заряд.
    
    Если отношения заряда / полной сумме зарядов <= 0.0001 - нашли уровень Ферми.
    """

    for i in range(1, 5):
        left = Efleft + step * (i - 1)
        right = Efleft + step * i

        n = calc_n(nc, left, right, 300.)
        # p = calc_p()
        print(n)
