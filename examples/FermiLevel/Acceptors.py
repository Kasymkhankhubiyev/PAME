from pame.ChargedParticlesInSemicondactor.AcceptorFermiLevel import *

if __name__ == '__main__':
    # поиск урвня Ферми для акцепторов
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=1e16)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=2 * 1e16)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=4 * 1e16)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=8 * 1e16)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=1e17)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=2 * 1e17)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=4 * 1e17)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=8 * 1e17)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=1e18)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=2 * 1e18)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=4 * 1e18)
    print(result)
    result = find_fermi_level(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=8 * 1e18)
    print(result)
    # result = calculate_charges(me=0.36, mh=0.36, t=300, Efpl=0.0, Efneg=0.57, Ec=1.12, Ev=0, Nd=1e19)
    # print(result)