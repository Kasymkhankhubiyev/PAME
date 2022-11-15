import matplotlib.pyplot as plt

def volt_amper_characteristics() -> None:

    def Id(ud: float, vp: float, ug: float):
        Ip = 0.1  # A ток отсечки
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
    plt.legend('best')
    plt.savefig('volt_amper_charact')


def g_m():
    pass


def g_d():
    pass