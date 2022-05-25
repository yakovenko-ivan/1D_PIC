import numpy as np

sor_eps = 1e-05
L = 10
cells = 127
nodes = cells + 1
dx = L / cells
rho = np.ones(nodes)
eps0 = 1.0


def SOR(rho):
    omega = 2.0 / (1 + np.pi / nodes)
    phi = np.array([0.0 for i in range(nodes)])
    rhs = -np.copy(rho) * dx**2 / eps0 
    for k in range (10000):
        phinew = np.copy(phi)        
        for i in range(nodes - 1):
            nxt = (i+1) if (i < nodes - 2) else 0
            prv = (i-1) if (i > 0) else (nodes - 2)
            phinew[i] = (1.0 - omega) * phi[i] + (omega / -2.0)*(rhs[i] - phinew[prv] - phi[nxt])
        if k % 25 == 0:
            if np.max(np.abs(phinew - phi)) < sor_eps:
                phi = phinew
                break        
        phi = np.copy(phinew)       
    phi[nodes - 1] = phi[0]  
    return phi

def field_on_nodes(phi):
    efield = np.array([0.0 for i in range(nodes)])
    for i in range(nodes):
        nxt = (i+1) if (i < nodes - 1) else 0
        prv = (i-1) if (i > 0) else (nodes - 1)
        efield[i] = (phi[prv] - phi[nxt]) / (2* dx)
    return efield

phi = SOR(rho)
e = field_on_nodes(phi)
print(e)
  