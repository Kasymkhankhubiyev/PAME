import matplotlib.pyplot as plt
# Ip = 0.1  # A ток отсечки

def volt_amper_characteristics(Ip: float) -> None:

    def Id(ud: float, vp: float, ug: float):
        res = Ip * (3*ud/vp - 2 * ((ud + ug)**1.5)/(vp**1.5) - (ug/vp)**1.5)
        return res

    fig = plt.Figure()
    for i in range(1, 10):
        Id_array, vp_array = [], []
        ud = i
        for j in range(1, 10):
            vp = ud*j/10
            id = Id(ud=ud, vp=vp, ug=1.)
            vp_array.append(vp)
            Id_array.append(id)
        plt.plot(vp_array, Id_array, label=f'ud={ud}')
    # plt.legend('best')
    plt.savefig('volt_amper_charact')


def g_m(Ip: float, ud: float, vp: float, ug: float):
    gm = -Ip * (3 * (ud + ug)**0.5/vp**1.5 - 3 * ug**0.5 / vp**1.5)
    return gm


def g_d(ud: float, vp: float, ug: float):
    gd = Ip * (3/vp - 3 * (ud + ug)**0.5 / vp**1.5)
    return gd
