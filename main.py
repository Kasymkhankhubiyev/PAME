from task1 import runfile
from ChargedParticlesInSemicondactor.DonorFermiLevel import *

if __name__ == '__main__':

    result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=10e17)
    print(result)

    # TODO сделать runfile, в который просто все передадим, а тот все соберет
    """
    Нужно просто разбить линейно область множества энергий Ферми
    передаем в расчет точку и считаем заряд.
    
    Если отношения заряда / полной сумме зарядов <= 0.0001 - нашли уровень Ферми.
    """
