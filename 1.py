import numpy as np
from matplotlib import pyplot as plt

wt = 1
lam = [0.2 + 0.05*i for i in range(1, 101)]
c = 1  # replication cost
dt = 5  # threshold
Pz = 0  # probability of completing the whole task
x2 = 0
Max = 0


def qiu_z(u, zd, zu, m, z):
    for i in range(100):
        M = -10000000
        n = m/c
        a = [1, n, zd[i], zu[i]]
        x1 = 0
        for j in range(len(a)):
            if a[j] < 1:
                continue
            if np.log(1 - np.exp(-a[j] * lam[i] * (dt - wt)))-u*c*a[j] > M:
                M = np.log(1 - np.exp(-a[j] * lam[i] * (dt - wt))) - u * c * a[j]
                x1 = j
        z[i] = a[x1]
    print(z)

def bm_sf_u(m):
    ur = 0.8
    ul = 0
    u = 0.35
    Epsilon = 0.00005
    d = 1
    z = 1 / (np.array(lam)*(dt - wt)) * np.log(1 + np.array(lam)*(dt - wt) / (u * c))
    zd = np.trunc(z).astype(int)
    zu = zd+1
    print("original:")
    qiu_z(u, zd, zu, m, z)
    while c*np.sum(z) != m and (ur - ul) >= Epsilon:
        if c*np.sum(z) > m:
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
        if c*np.sum(z) == m:
            print("exit the loop")
            break
    print("the best u:%f" % u)
    print("the best allocation method:")
    print(z)


def dw_f_la(m):
    global Pz, x2, Max
    x2 = 0
    Max = 0
    R = np.ones(100)
    logP = np.log(1 - np.exp(-R * lam * (dt - wt)))
    for j in range(1, m-100+1):
        for i in range(100):
            logPi1 = np.log(1 - np.exp(-(R[i]+1) * lam[i] * (dt - wt)))
            dif = logPi1-logP[i]
            if dif > Max:
                x2 = i
                Max = dif
        Max = 0
        logP[x2] = np.log(1 - np.exp(-(R[x2] + 1) * lam[x2] * (dt - wt)))
        R[x2] += 1
    P = 1 - np.exp(-R * np.array(lam) * (dt - wt))
    Pz = 1-np.prod(P)
    return R

if __name__ == '__main__':
    bm_sf_u(150)
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
        ax[a].set_title("Fig.a%d. Best distribution(M=%d)" % (a+1,110+a*20), fontsize=15)
        ax[a].bar(x, height=dw_f_la(110+a*20), width=1, color='red', edgecolor='black', align='edge')
    ax4.set_xlabel("M", fontsize=12)
    ax4.set_ylabel("Task outage probability", fontsize=15)
    ax4.set_title("Fig.a4. Cost and delay trade-off", fontsize=15)
    y4 = np.zeros(101)
    for a in range(100, 201):
        dw_f_la(a)
        y4[a-100] = Pz
    ax4.plot(np.arange(100, 201), y4, "r .")
    plt.show()
