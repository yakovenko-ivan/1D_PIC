import numpy as np
from math import *
import matplotlib.pyplot as plt

ng, epsi, l, T, nth, A1, A2, iw, e0, w0 = 1024, 4*np.pi, 1, 1, 10, 2, 1, 2, 1, 1
dx = l/(ng-1)
dt = T/(nth-1)
x = np.linspace(0, l+dx, ng+1)
global ng, epsi, l, T, dx, dt, nth, A1, A2, iw, e0, w0

class Fields:
  def __init__(self):
    self.ese = np.zeros(nth)
    self.e = np.zeros(ng+1)
    self.rho = np.zeros(ng+1)
    self.rho0 = np.zeros(ng+1)
    #self.rho0[ng//2] = 1
    self.rho0 = np.cos(2*np.pi*x/(0.2*l))
    for k in range(ng+1):
        self.rho[k] = self.rho0[k]
    self.phi = np.zeros(ng+1)
    self.ael = 1
    
    
  def compute_fields(self, ith):

    # ratio phik/rhok
    ng2 = ng // 2
    sm = []
    ksqi = []
    for k in range(ng2+1):
        kdx = (2.0*pi / ng) * (k + 1)
        sm.append(exp(A1*sin(kdx)**2 - A2*tan(kdx)**4))
        #sm.append(1)
        ksqi.append(epsi / (2.0*sin(kdx/2.0)/dx)**2 * sm[k]**2)
    
    # transform charge density
    self.rho[0] = self.rho[0] + self.rho[ng]
    self.rho[ng] = self.rho[0]
    hdx = 0.5 * dx
    rhok = self.rho*hdx
    rhok = np.fft.rfft(rhok)
    rhok[0] = 0.0
    
    # calculate phik and field energy
    eses = 0.0
    phik = np.zeros(ng, dtype=np.complex_)
    phik = rhok*ksqi
    for k in range(1,ng2):
        eses = eses + rhok[k] * phik[k]
       
    phik[ng2] = ksqi[ng2-1] * rhok[ng2]
    #self.ese[ith + 1] = (2.0 * eses + rhok[ng2] * phik[ng2]) / (2.0 * l)
    self.ese[ith] = (2.0 * eses + rhok[ng2] * phik[ng2]) / (2.0 * l)
    rhok[ng2] = sm[ng2-1] * rhok[ng2]
    
    # inverse transform
    li = 1.0/l
    for k in range(ng2+1):
        rhok[k] = rhok[k] * li
        phik[k] = phik[k] * li
    self.rho = np.fft.irfft(rhok).real
    self.phi = np.fft.irfft(phik).real
    self.rho = np.append(self.rho, self.rho[0])
    self.phi = np.append(self.phi, self.phi[0])
    
    # uniform field
    time = ith*dt
    e0t = e0 * np.cos(w0*time)
    
    # centred difference
    if (iw == 1 or iw == 2):
        hdxi = 0.5/dx
        for j in range(1, ng):
            self.e[j] = (self.phi[j-1] - self.phi[j+1]) * hdxi + e0t
        self.e[0] = (self.phi[ng] - self.phi[1]) * hdxi + e0t
        self.e[ng] = self.e[0]
        
    else:
        dxi = 1.0/dx
        for j in range(ng):
            self.e[j] = (self.phi[j] - self.phi[j+1]) * dxi + e0t
        self.e[0] = (self.phi[ng] - self.phi[1]) * hdxi + e0t
        self.e[ng] = self.e[0]
        
    # clear out old charge density
    for k in range(ng + 1):
        self.rho[k] = self.rho0[k]
    self.rho[ng] = 0.0
    
    # electric field has not been renormalized yet
    self.ael = 1


field = Fields()

for ith in range(nth):
    field.compute_fields_rfft(ith)
    
#print(field.e)
#plt.plot(x, field.e)
plt.plot(x, field.phi/max(abs(field.phi)))
plt.plot(x, field.rho/max(abs(field.rho)))