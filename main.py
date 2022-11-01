from task1 import runfile
from ChargedParticlesInSemicondactor.DonorFermiLevel import *
from SemiconCurrent import CurrentCalculus

if __name__ == '__main__':

    # result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.57, Efneg=1.12, Ec=1.12, Ev=0, Nd=10e17)
    # print(result)

    Nd = 5e17
    Na = 1e18

    js = CurrentCalculus.count_Js(t=300., Nd=Nd, Na=Na)
    print(js)



