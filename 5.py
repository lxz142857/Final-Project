import numpy as np
from matplotlib import pyplot as plt
import math

wt = 1
lam = [0.2 + 0.1*i for i in range(1, 51)]
xm = 1
alp = [1.2 + 0.1*i for i in range(1, 51)]
c = 1
dt = 5
Exr = wt + 1/np.array(lam)  # expected completion time of SED tasks in descending order
Exb = xm + xm/(np.array(alp)-1)  # expected completion time of Pareto tasks in descending order
color = np.ones(100).astype(int)  # color = 1--red--SED color=0--blue--Pareto
ind = np.zeros(100).astype(int)  # index of subtasks in decreasing order
Ex = np.zeros(100)
Pz = 0
x2 = 0
Max = 0


def i_d_o():
    r1 = 0
    b1 = 0
    for i in range(100):
        if r1 == 50 and b1 < 50:
            s = 50 - b1
            for j in range(s):
                color[50+b1] = 0
                Ex[i] = Exb[b1]
                b1 += 1
            break
        if b1 == 50 and r1 < 50:
            s = 50 - r1
            for j in range(s):
                color[50 + r1] = 1
                Ex[i] = Exr[r1]
                r1 += 1
            break
        if Exr[r1] >= Exb[b1]:
            color[i] = 1
            Ex[i] = Exr[r1]
            r1 += 1
        else:
            color[i] = 0
            Ex[i] = Exb[b1]
            b1 += 1
    r2 = 0
    b2 = 0
    for i in range(100):
        if color[i] == 0:
            ind[i] = b2
            b2 += 1
        else:
            ind[i] = r2
            r2 += 1


def qiu_z(u, zrd, zru, zbd, zbu, m, zr, zb):
    for i in range(50):
        M = -10000000
        n = m/c
        a = [1, n, zrd[i], zru[i]]
        x1 = 0
        for j in range(len(a)):
            if a[j] < 1:
                continue
            if np.log(1 - np.exp(-a[j] * lam[i] * (dt - wt)))-u*c*a[j] > M:
                M = np.log(1 - np.exp(-a[j] * lam[i] * (dt - wt))) - u * c * a[j]
                x1 = j
        zr[i] = a[x1]
    for j in range(50):
        M = -10000000
        n = m/c
        a = [1, n, zbd[j], zbu[j]]
        x1 = 0
        for i in range(len(a)):
            if a[i] < 1:
                continue
            if np.log(1 - math.pow(xm/dt, a[i] * alp[j]))-u*c*a[i] > M:
                M = np.log(1 - math.pow(xm/dt, a[i] * alp[j]))-u*c*a[i]
                x1 = i
        zb[j] = a[x1]


def bm_sf_u(m):
    ur = 0.2
    ul = 0
    u = 0.01
    Epsilon = 0.00005
    d = 1
    zr = 1 / (np.array(lam)*(dt - wt)) * np.log(1 + np.array(lam)*(dt - wt) / (u * c))
    zb = np.log(u * c / (u * c - np.array(alp) * math.log(xm / dt))) / (np.array(alp) * math.log(xm / dt))
    zrd = np.trunc(zr).astype(int)
    zru = zrd+1
    zbd = np.trunc(zb).astype(int)
    zbu = zbd+1
    qiu_z(u, zrd, zru, zbd, zbu, m, zr, zb)
    z = np.append(zr, zb)
    while c*np.sum(z) != m and (ur - ul) >= Epsilon:
        if c*np.sum(z) > m:
            ul = u
        else:
            ur = u
        u = (ur+ul)/2
        zr = 1 / (np.array(lam) * (dt - wt)) * np.log(1 + np.array(lam) * (dt - wt) / (u * c))
        zb = np.log(u * c / (u * c - np.array(alp) * math.log(xm / dt))) / (np.array(alp) * math.log(xm / dt))
        zrd = np.trunc(zr).astype(int)
        zru = zrd + 1
        zbd = np.trunc(zb).astype(int)
        zbu = zbd + 1
        print("search u for the %d time(s)" % d)
        d += 1
        qiu_z(u, zrd, zru, zbd, zbu, m, zr, zb)
        z = np.append(zr, zb)
        print(z)
        if c*np.sum(z) == m:
            print("exit the loop")
            break
    for i in range(100):
        if color[i] == 1:
            z[i] = zr[ind[i]]
        else:
            z[i] = zb[ind[i]]
    print("the best u:%f" % u)
    print("the best allocation method：")
    print(z)


def dw_f_la(m):
    global Pz, x2, Max
    x2 = 0
    Max = 0
    Rr = np.ones(50)
    Rb = np.ones(50)
    R = np.ones(100)
    logPr = np.log(1 - np.exp(-Rr * lam * (dt - wt)))
    logPb = np.log(1 - np.power(xm / dt, Rb * np.array(alp)))
    for j in range(1, m-100+1):
        for i in range(100):
            if color[i] == 1:
                logPr1 = np.log(1 - np.exp(-(Rr[ind[i]] + 1) * lam[ind[i]] * (dt - wt)))
                dif = logPr1 - logPr[ind[i]]
            else:
                logPb1 = np.log(1 - np.power(xm / dt, (Rb[ind[i]] + 1) * alp[ind[i]]))
                dif = logPb1 - logPb[ind[i]]
            if dif > Max:
                x2 = i
                Max = dif
        Max = 0
        if color[x2] == 1:
            logPr[ind[x2]] = np.log(1 - np.exp(-(Rr[ind[x2]] + 1) * lam[ind[x2]] * (dt - wt)))
            Rr[ind[x2]] += 1
        else:
            logPb[ind[x2]] = np.log(1 - np.power(xm / dt, (Rb[ind[x2]] + 1) * alp[ind[x2]]))
            Rb[ind[x2]] += 1
        R[x2] += 1
    P1 = 1 - np.exp(-Rr * np.array(lam) * (dt - wt))
    P2 = 1 - np.power(xm / dt, Rb * np.array(alp))
    P = np.append(P1, P2)
    Pz = 1-np.prod(P)
    return R


if __name__ == '__main__':
    i_d_o()
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
        ax[a].set_ylim(0, 5)
        ax[a].set_xticks(np.linspace(0, 100, 11))
        ax[a].set_yticks(np.linspace(1, 5, 5))
    ax4.set_xlim(100, 200)
    ax4.set_ylim(0, 1)
    ax4.set_xticks(np.linspace(100, 200, 11))
    ax4.set_yticks(np.linspace(0, 1, 11))
    x = np.arange(0, 100, 1)
    for a in range(3):
        ax[a].set_xlabel("Subtask index i", fontsize=12)
        ax[a].set_ylabel("The number of replications for subtask i", fontsize=12)
        ax[a].set_title("Fig.e%d. Best distribution(M=%d)" % (a+1,110+a*20), fontsize=15)
        ax[a].bar(x, height=dw_f_la(110+a*20), width=1,
                  color=np.where(color[x] == 1, 'r', 'b'), edgecolor='black', align='edge')
    ax4.set_xlabel("M", fontsize=12)
    ax4.set_ylabel("Task outage probability", fontsize=15)
    ax4.set_title("Fig.e4. Cost and delay trade-off", fontsize=15)
    y4 = np.zeros(101)
    for a in range(100, 201):
        dw_f_la(a)
        y4[a-100] = Pz
    ax4.plot(np.arange(100, 201), y4, "r .")
    plt.show()
