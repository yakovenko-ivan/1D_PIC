il = 1
iu = 1


# тут не особо понятно по какому i идет итерирование
def f100():
    v1s = 0
    v2s = 0
    for i in range(il, iu):
        f101(i)
    j = x[i] + 0.5
    vo = vs[i]
    vn = vo + a[j + 1]
    v1s = v1s + vn
    v2s = v2s + vn * vo


def f101(i):
    vx[i] = vn
    p = p + m * v1s * dxdt
    ke = ke + 0.5 * m * v2s * dxdt ** 2


def f200():
    if t != 0:
        f250()
    v1s = 0
    v2s = 0
    if vec:
        f2000()


def f2009(i):
    for i in range(il, iu):
        f201(i)
    j = x[i]
    vo = vx[i]
    vn = no + a[j + 1] + (x[i]- j) * (a[j + 2] - a[j + 1])
    v1s = v1s + vn
    v2s = v2s + vo * vn


def f201(i):
    vx[i] = vn
    p = p + m * v1s * dxdt
    ke = ke + 0.5 * m * v2s * dxdt ** 2


def f300():
    v1s = 0
    v2s = 0
    for i in range(il, iu):
        f301(i)
    j = x[i]
    vo = vx[i]
    vn = vo + a[j + 1]
    v1s = v1s + vn
    v2s = v2s + vn * vo


def f301(i):
    vx[i] = vn
    p = p + m * v1s * dxdt
    ke = ke + 0.5 * m * v2s * dxdt ** 2


def f250():
    s = 2 * t / (1 + t ** 2)
    v2s = 0
    if vec:
        f2500()


def f2509(i):
    for i in range(il, iu):
        f251(i)
    j = x[i]
    aa = a[j+1] + (x[i] - j) * (a[j+2] - a[j+1])
    vyy = vy[i]
    vxx = vx[i] - t * vyy + aa
    vyy = vyy + s * vxx
    vxx = vxx - t * vyy
    v2s = v2s + vxx * vxx * vyy * vyy
    vx[i] = vxx + aa


def f251(i):
    vy[i] = vyy
    ke = ke + 0.5 * m * v2s * dxdt ** 2


def f2000(i):
    for i in range(1, 64):
        f2004(i)
    v1si[i] = 0


def f2004():
    v2si[i] = 0
    for j in range(il, iu - 63, 64):
        f2006(j)
    for i in range(0, 63):
        f2001(i)


def f2001(i):
    ji[i + 1] = x[i + j]
    for i in range(1, 64, 2):
        f2002(i)
    al[i + 1] = a[ji[i + 1] + 1]
    ar[i + 1] = a[ji[i + 1] + 2]
    al[i] = a[ji[i] + 1]


def f2002(i):
    ar[i] = a[ji[i] + 2]
    for i in range(0, 63):
        vni[i + 1] = vx[i + j] + al[i + 1] + (x[i + j] - ji[i + 1]) * (ar[i + 1] - al[i + 1])
        v1si[i + 1] = v1si[i + 1] + vni[i + 1]
        v2si[i + 1] = v2si[i + 1] + vx[i + j] * vni[i + 1]


def f2005():
    vx[i + j] = vni[i + 1]


def f2006():
    il = il + 64
    for i in range(1, 64):
        f2007(i)
    v1s = v1s + v1si[i]


def f2007(i):
    v2s = v2s + v2si[i]
    f2009(i)


def f2500():
    for i in range(1, 64):
        f2504(i)


def f2504(i):
    v2si[i] = 0
    for j in range(il, iu - 63, 64):
        f2506(j)
    for i in range(0, 63):
        f2501(i)


def f2501(i):
    ji[i + 1] = x[i + j]
    for i in range(1, 64, 2):
        f2502()
    al[i + 1] = a[ji[i + 1] + 1]
    ar[i + 1] = a[ji[i + 1] + 2]
    al[i] = a[ji[i] + 1]


def f2502():
    ar[i] = a[ji[i] + 2]
    for i in range(0, 63):
        f2505()
    aai[i + 1] = al[i + 1] + (x[i+j] - ji[i + 1]) * (ar[i+1] - al[i+1])
    vx[i + j] = vx[i + j] - t * vy[i + j] + aai[i + 1]
    vy[i + j] = vy[i + j] + s * vx[i + j]
    vx[i + j] = vx[i + j] - t * vy[i + j]
    v2si[i + j] = v2si[i + 1] + vx[i + j] ** 2 + vy[i + j] ** 2


def f2505():
    vx[i + j] = vx[i + j] + aai[i + 1]


def f2506(j):
    il = il + 64
    for i in range(1, 64):
        f2507(i)


def f2507(i):
    v2s = v2s + v2si[i]
    f2509(i)




