from task1 import runfile
from ChargedParticlesInSemicondactor.DonorFermiLevel import *
from SemiconCurrent import CurrentCalculus

if __name__ == '__main__':

    # поиск урвня Ферми
    # result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=10e17)
    # print(result)

    # рассчет тока в полупроводнике в pn-переходе
    # Nd = 5e17
    # Na = 1e18
    #
    # js = CurrentCalculus.count_Js(t=300., Nd=Nd, Na=Na)
    # print(js)

    # рассчет уровня энергии в периодической яме
    runfile.run()



