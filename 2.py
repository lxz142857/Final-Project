import numpy as np
from matplotlib import pyplot as plt

wt = 1
lam = [0.2 + 0.05*i for i in range(1, 101)]
N = 100
rho = N/np.sum((wt+1/np.array(lam)))  # coefficient of replication cost
c = rho*(wt+1/np.array(lam))
dt = 5
Pz = 0


def qiu_z(u, zd, zu, m, z):
    for i in range(100):
        M = -10000000
        n = np.trunc(m/c).astype(int)
        a = [1, n[i], zd[i], zu[i]]
        x1 = 0
        for j in range(len(a)):
            if 1 - np.exp(-a[j] * lam[i] * (dt - wt)) > 0:
                if np.log(1 - np.exp(-a[j] * lam[i] * (dt - wt)))-u*c[i]*a[j] > M:
                    M = np.log(1 - np.exp(-a[j] * lam[i] * (dt - wt))) - u * c[i] * a[j]
                    x1 = j
        z[i] = a[x1]
    print(z)


def bm_sf_u(m):
    global Pz
    ur = 0.3
    ul = 0
    u = 0.000001
    Epsilon = 0.000000001
    d = 1
    z = 1 / (np.array(lam)*(dt - wt)) * np.log(1 + np.array(lam)*(dt - wt) / (u * c))
    zd = np.trunc(z).astype(int)
    zu = zd+1
    qiu_z(u, zd, zu, m, z)
    while np.sum(c*z) != m and (ur - ul) >= Epsilon:
        if np.sum(c*z) > m:
            ul = u
        else:
            ur = u
        u = (ur+ul)/2
        z = 1 / (np.array(lam)*(dt - wt)) * np.log(1 + np.array(lam)*(dt - wt) / (u * c))
        zd = np.trunc(z).astype(int)
        zu = zd + 1
        print("search u for the %d time(s)" % d)
        d += 1
        qiu_z(u, zd, zu, m, z)
        if np.sum(c*z) == m:
            print("exit the loop")
            break
    print("the best u:%f" % u)
    print("the best allocation method:")
    print(z)
    P = 1 - np.exp(-z * np.array(lam) * (dt - wt))
    Pz = 1 - np.prod(P)
    return z


if __name__ == '__main__':
    fig = plt.figure(num=1, figsize=(15, 12))
    plt.rcParams["font.family"] = "Microsoft Yahei"
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    ax = [ax1, ax2, ax3, ax4]
    for a in range(3):
        ax[a].set_xlim(0, 100)
        ax[a].set_ylim(0, 6)
        ax[a].set_xticks(np.linspace(0, 100, 11))
        ax[a].set_yticks(np.linspace(1, 6, 6))
    ax4.set_xlim(100, 200)
    ax4.set_ylim(0, 1)
    ax4.set_xticks(np.linspace(100, 200, 11))
    ax4.set_yticks(np.linspace(0, 1, 11))
    x = np.arange(0, 100, 1)
    for a in range(3):
        ax[a].set_xlabel("Subtask index i", fontsize=12)
        ax[a].set_ylabel("The number of replications for subtask i", fontsize=12)
        ax[a].set_title("Fig.b%d. Best distribution(M=%d)" % (a+1,110+a*20), fontsize=15)
        ax[a].bar(x, height=bm_sf_u(110+a*20), width=1, color='red', edgecolor='black', align='edge')
    ax4.set_xlabel("M", fontsize=12)
    ax4.set_ylabel("Task outage probability", fontsize=15)
    ax4.set_title("Fig.b4. Cost and delay trade-off", fontsize=15)
    y4 = np.zeros(101)
    for a in range(100, 201):
        bm_sf_u(a)
        y4[a-100] = Pz
    ax4.plot(np.arange(100, 201), y4, "r .")
    plt.show()
