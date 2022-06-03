#!/usr/bin/env python
# coding: utf-8

# In[11]:


import numpy as np


# In[25]:


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


# In[26]:


def SETRHO(il, iu, q, rhos):
    
    qdx = q/dx
    if(il == 1):
        for j in range(1, ng):
            rho[i] = rho0
        rho[ng+1] = 0
        dxi = 1.0/dx
        xn = ng
    
    rho0 = rho0 - rhos
    for j in range(1, ng):
        rho[j] = rho[j] - rhos
    
    if iw == 1:
        for i in range(il, iu):
            x[i] = x[i] * dxi
            if(x[i] < 0): x[i] += xn
            if(x[i] >=0): x[i] -= xn
            j = x[i] + 0.5
            rho[j+1] += qdx    
    if iw == 2 or iw == 3:
        for i in range(il, iu):
            x[i] = x[i] * dxi
            if(x[i] < 0): x[i] += xn
            if(x[i] >=0): x[i] -= xn
            j = x[i]
            drho = qdx*(x[i] - j)
            rho[j+1] += - drho + qdx
            rho[j+2] += drho


# In[27]:


def SETV(il, iu, q, m, t, p):
    
    dtdx = dt/dx
    if(t != 0):
        c = 1.0/np.sqrt(1.0 + t*t)
        s = c*t
        for i in range(il, iu):
            vx[i] = c*vxx + s*vy[i]
            vy[i] = -s*vxx + c*vy[i]
            vy[i] =  vy[i]*dtdx
    for i in range(il, iu):
        vx[i] = vx[i]*dtdx
    
    # eleectric impulse to go back 1/2 time step
    data = []
    # Вызываем процедуру, написанную второй группой
    accel(ilp, iup, -0.5*q, m, 0.0, p, data)

