import numpy as np


ilp = 10
iup = 10 * 64 * 10
q = 1
m = 1
t = 1
p = 1
ke = 1
dx, dt = 1, 1
vn = 1
vec = 1

dxdt = dx / dt
ae = (q/m) * (dt / dxdt)


a = np.empty(shape=64, dtype=np.float64)
ji = np.empty(shape=64, dtype=np.float64)
al = np.empty(shape=64, dtype=np.float64)
ar = np.empty(shape=64, dtype=np.float64)
vni = np.empty(shape=64, dtype=np.float64)
v1si = np.empty(shape=64, dtype=np.float64)
v2si = np.empty(shape=64, dtype=np.float64)
aai = np.empty(shape=64, dtype=np.float64)
vx = []
vy = []
x = []
il = ilp
iu = iup


def f100():
    v1s = 0
    v2s = 0
    for i in range(il, iu):
        vx[i] = vn
        p = p + m * v1s * dxdt
        ke = ke + 0.5 * m * v2s * dxdt ** 2
    j = x[i] + 0.5
    vo = vx[i]
    vn = vo + a[j + 1]
    v1s = v1s + vn
    v2s = v2s + vn * vo


def f200():
    if t != 0:
        s = 2 * t / (1 + t ** 2)
        v2s = 0
        if vec:
            for i in range(1, 64):
                v2si[i] = 0
                for j in range(il, iu - 63, 64):
                    il = il + 64
                    for i in range(1, 64):
                        v2s = v2s + v2si[i]
                        for i in range(il, iu):
                            vy[i] = vyy
                            ke = ke + 0.5 * m * v2s * dxdt ** 2
                        j = x[i]
                        aa = a[j + 1] + (x[i] - j) * (a[j + 2] - a[j + 1])
                        vyy = vy[i]
                        vxx = vx[i] - t * vyy + aa
                        vyy = vyy + s * vxx
                        vxx = vxx - t * vyy
                        v2s = v2s + vxx * vxx * vyy * vyy
                        vx[i] = vxx + aa
                for i in range(0, 63):
                    ji[i + 1] = x[i + j]
                    for i in range(1, 64, 2):
                        ar[i] = a[ji[i] + 2]
                        for i in range(0, 63):
                            vx[i + j] = vx[i + j] + aai[i + 1]
                        aai[i + 1] = al[i + 1] + (x[i + j] - ji[i + 1]) * (ar[i + 1] - al[i + 1])
                        vx[i + j] = vx[i + j] - t * vy[i + j] + aai[i + 1]
                        vy[i + j] = vy[i + j] + s * vx[i + j]
                        vx[i + j] = vx[i + j] - t * vy[i + j]
                        v2si[i + j] = v2si[i + 1] + vx[i + j] ** 2 + vy[i + j] ** 2
                    al[i + 1] = a[ji[i + 1] + 1]
                    ar[i + 1] = a[ji[i + 1] + 2]
                    al[i] = a[ji[i] + 1]
    v1s = 0
    v2s = 0
    if vec:
        for i in range(1, 64):
            v2si[i] = 0
            for j in range(il, iu - 63, 64):
                il = il + 64
                for i in range(1, 64):
                    v2s = v2s + v2si[i]
                    for i in range(il, iu):
                        vx[i] = vn
                        p = p + m * v1s * dxdt
                        ke = ke + 0.5 * m * v2s * dxdt ** 2
                    j = x[i]
                    vo = vx[i]
                    vn = vo + a[j + 1] + (x[i] - j) * (a[j + 2] - a[j + 1])
                    v1s = v1s + vn
                    v2s = v2s + vo * vn
                v1s = v1s + v1si[i]
            for i in range(0, 63):
                ji[i + 1] = x[i + j]
                for i in range(1, 64, 2):
                    ar[i] = a[ji[i] + 2]
                    for i in range(0, 63):
                        vni[i + 1] = vx[i + j] + al[i + 1] + (x[i + j] - ji[i + 1]) * (ar[i + 1] - al[i + 1])
                        v1si[i + 1] = v1si[i + 1] + vni[i + 1]
                        v2si[i + 1] = v2si[i + 1] + vx[i + j] * vni[i + 1]
                al[i + 1] = a[ji[i + 1] + 1]
                ar[i + 1] = a[ji[i + 1] + 2]
                al[i] = a[ji[i] + 1]
        v1si[i] = 0


def f300():
    v1s = 0
    v2s = 0
    for i in range(il, iu):
        vx[i] = vn
        p = p + m * v1s * dxdt
        ke = ke + 0.5 * m * v2s * dxdt ** 2
    j = x[i]
    vo = vx[i]
    vn = vo + a[j + 1]
    v1s = v1s + vn
    v2s = v2s + vn * vo