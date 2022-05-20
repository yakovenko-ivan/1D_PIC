def init(ins: int, ins_1: int, ms: float, qs: float, ts:float, nms: int, rho0:float):
    a = 1
    return a


def setv(ins: int, ins_1: int, qs: float, ms: float, ts: float, pxs: float):
    a = 1
    return a


def fields(var):
    a = 1
    return a


def plotxv(var_1, ins: int, vl: float, vu: float):
    a = 1
    return a


def plotfvx(var_1, ins: int, vl: flaot, vu: float, qs: float, var_2):
    a = 1
    return a


def histry():
    a = 1
    return a


def last():
    a = 1
    return a


if __name__ == '__main__':
    print('This is ES1 the program used for one-dimensional electrostatic problems solving')
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--length', help='physical length of system', type=float, default=6.283185307177958)
    parser.add_argument('-nsp', '--number_sp', help='number of particle species', type=int, default=1)
    parser.add_argument('-dt', '--delta_t', help='time step', type=float, default=0.2)
    parser.add_argument('-nt', '--number_t', help='number of time steps to run (ending time = nt * dt', type=int,
                        default=150)
    parser.add_argument('-ng', '--number_g', help='number of spacial cells, must be power of 2', type=int, default=32)
    parser.add_argument('-iw', '--mover', help='mover algorithm selector, see accel -> [100, 200, 300] and move',
                        type=int, default=2)
    parser.add_argument('-vec', '--vectorization', help='to select vectorized option where possible', type=bool,
                        default=True)
    parser.add_argument('-epsi', '--epsilon', help='multiplier in poisson equation, 1 for rationalized units',
                        type=float, default=1.0)
    parser.add_argument('-a1', '--boost_1', help='filed smoothing boost', type=float, default=0.0)
    parser.add_argument('-a2', '--boost_2', help='mid-range boost', type=float, default=0.0)
    parser.add_argument('-е0', '--field_amplitude', help='add uniform electric field e0*cos(w0*time)', type=float,
                        default=0.0)
    parser.add_argument('-w0', '--field_frequency', help='add uniform electric field e0*cos(w0*time)', type=float,
                        default=0.0)
    parser.add_argument('-irho', '--plot_interval', help='plotting interval for rho (charge density)', type=int)
    parser.add_argument('-irhos', '--plot_interval_smooth', help='plotting interval for smoothed density', type=int)
    parser.add_argument('-iphi', '--plot_interval_potential', help='plotting interval for phi (potential)', type=int)
    parser.add_argument('-ie', '--plot_interval_field', help='plotting interval for electric field', type=int)
    parser.add_argument('-ixvx', '--plot_interval_x_vx', help='plotting interval for x vs vx phase space', type=int)
    parser.add_argument('-ifvx', '--plot_interval_f_vx', help='plotting interval for f(vx) distribution'
                                                              '> 0 gives linear, <0 gives semi-log',
                        type=int, default=0)
    parser.add_argument('-ixvy', '--plot_interval_x_vy', help='plotting interval for x vs vx phase space', type=int)
    parser.add_argument('-mplot', '--m_plots', help='fourier mode numbers to plot', type=int, default=0)

    args = parser.parse_args()
    l = args.length
    nsp = args.number.sp
    dt = args.delta_t
    nt = args.numer_t
    ng = args.number_g
    iw = args.mover
    vec = args.vectorization
    epsi = args.epsilon
    a1 = args.boost_1
    a2 = args.boost_2
    e0 = args.field_amplitude
    w0 = args.field_frequency
    irho = args.plot_interval
    irhos = args.plot_interval_smooth
    iphi = args.plot_interval_potential
    ie = args.plot_interval_field
    ixvx = args.plot_interval_x_vx
    ifvx = args.plot_interval_f_vx
    ixvy = args.plot_interval_x_vy
    mplot = args.m_plots



#   from cliche mfield
    ngmax = 256
    ng1m = ngmax + 1
    rho = np.empty(shape=ng1m)
    phi = np.empty(shape=ng1m)
    e = np.empty(shape=ng1m)

#   from cliche mptcl
    x = np.empty(shape=8192)
    vx = np.empty(shape=8192)
    yx = np.empty(shape=8192)

#   from cliche mtime
    nth = 500
    mmax = 10
    nspm = 3
    nth1, nth2, nspm1 = nth + 1, nnth + 2, nspm + 1
    ese = np.empty(shape=nth1)
    kes = np.empty(shape=(nth1, nspm), dtype=float)
    pxs = np.empty(shape=(nth2, nspm), dtype=float)
    nms = np.empty(shape=nspm)
    # mplot = mmax
    if mplot == 0:
        mplot = mmax
    esem = np.empty(shape=(nth1, mmax))
    ms = np.empty(shape=nspm, dtype=float)
    qs = np.empty(shape=nspm, dtype=float)
    ts = np.empty(shape=nspm, dtype=float)
    ntp = 100
    rho0 = 0
    it, time, ith, ithl = 0, 0, 0, 0
    ins = [1]


    def accel(ilp: int, iup: int, q: float, m: float, t: float, p: float, ke: float) -> (float, float):
        import numpy as np
        # This function is called to convert E to (q/m)E*dt^2/dx, called A, a computer variable. Advanced velocity
        # one time step, using weighted E (or A). Calculate momentum and kinetic energy. Being repeated for each
        # species.
        a = np.empty(shape=64, dtype=np.float64)
        v1si = np.empty(shape=64, dtype=np.float64)
        v2si = np.empty(shape=64, dtype=np.float64)
        # vx = np.ndarray(shape=64, dtype=np.float64)
        il = ilp
        iu = iup

        dxdt = dx / dt
        ae = (q / m) * dt / dxdt
        if t != 0:
            ae = ae / 2
        if ae != ael:
            ng1 = ng + 1
            tem = ae / ael
            for j in range(1, ng1):
                a[j] = a[j] * tem
        if iw == 100:
            v1s = 0
            v2s = 0
            for i in range(il, iu):
                j = int(x[i] + 0.5)
                # print(f'{j}, {type(j)}')
                vo = vx[i]
                vn = vo + a[j + 1]
                v1s = v1s + vn
                v2s = v2s + vn * vo
                vx[i] = vn
            p = p + m * v1s * dxdt
            ke = ke + 0.5 * m * v2s * dxdt ** 2
            return p, ke
        elif iw == 200:
            if t == 0:
                v1s = 0
                v2s = 0
                if vec:
                    for i in range(1, 64):
                        v1si[i] = 0
                        v2si[i] = 0
            else:
                s = 2 * t / (1 + t ** 2)
                v2s = 0
                if vec:
                    for i in range(1, 64):
                        v2si[i] = 0
            for i in range(il, iu):
                j = int(x[i])
                vo = vx[i]
                vn = vo + a[j + 1] + (x[i] - j) * (a[j + 2] - a[j + 1])
                v1s = v1s + vn
                v2s = v2s * vo * vn
                vx[i] = vn
            p = p + m * v1s * dxdt
            ke = ke + 0.5 * m * v2s * dxdt ** 2
            return p, ke
        if iw == 300:
            v1s = 0
            v2s = 0
            for i in range(il, iu):
                j = int(x[i])
                vo = vx[i]
                vn = vo + a[j + 1]
                v1s = v1s + vn
                v2s = v2s + vn * vo
                vx[i] = vn
            p = p + m * v1s * dxdt
            ke = ke + 0.5 * m * v2s * dxdt ** 2
            print(f'{p}\n{ke}')
            return p, ke


    # here should be called histry
    dx = l/ng
    for var_is in range(1, nsp):
        # here should be called function init
        init(ins[var_is], ins[var_is + 1], ms[var_is], qs[var_is], ts[var_is], nms[var_is], rho0)

    fields(0)

    for var_is in range(1, nsp):
        # here should be called function setv
        setv(ins[var_is], ins[var_is + 1] - 1,qs[var_is],  ms[var_is], ts[var_is], pxs[1, var_is])

    def function_100(nsp, vmu):
        vl = 0.0
        vu = 0.0
        plotxv(1, ins[nsp+1] - 1, vl, vu)
        plotfvx(1, ins[2] - 1, vl, vu, qs[1], 32)
        if ts[1] != 0:
            pltvxy(1, ins[2] - 1, vmu)
        # p = 0
        # ke = 0
        for _is in range(1, nsp):
            p, ke = accle(ins[_is], ins[_is + 1] - 1, qs[_is], ms[_is], ts[_is], pxs[ith + 2,  _is], kes[ith + 1, _is])
            p = p + pxs[ith + 2, _is]
            ke = ke + kes[ith + 1, _is]
        for _ in range(1, nsp):
            move(ins[_], ins[_+1] - 1, qs[_])  # должна вернуть значение it
            te = ke + ese[ith+1]
        return

    if it < nt:
        if ith == nth:
            histry()
        it = it + 1
        time = it * dt
        ith = it - ithl
        fields(ith)
        function_100(nsp, vmu)
    histry()
    last()
