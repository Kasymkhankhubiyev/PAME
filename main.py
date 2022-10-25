from task1 import runfile
from ChargedParticlesInSemicondactor.chargedparticles import *

if __name__ == '__main__':

    # runfile.run()
    nc = calculate_charges(me=0.36, t=300, mh=0.81)
    print(nc)

    pl = (0.57 + 1.12)/2.

    n = calc_n(nc, pl, 1.12, 300)
    print(n)
