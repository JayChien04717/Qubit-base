import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
from qutip import *


#in MHz
g = 1.0 * 2 * pi # coupling strength
g1 = 0.5    # relaxation rate
g2 = 0.001       # dephasing rate
n_th = 10e-3       # bath temperature
f_q = 6

c_ops=[]
c_ops.append(sqrt(g1 * (1+n_th)) * sigmam())
c_ops.append(sqrt(g1 * n_th) * sigmam().dag())
c_ops.append(sqrt(g2) * sigmaz())


def const(x: np.array, lv: float):
    return lv * np.ones(x.size)


t1 = np.linspace(0,20,100)
t2 = np.linspace(21, 30, 100)
tlist = np.concatenate((t1,t2))


pulse = np.concatenate((const(t1,3), const(t2,0)))
# Hamiltonian
H = QobjEvo([0.5*sigmaz(), [sigmax(),pulse]], tlist=tlist)
H2 = QobjEvo([[0.5*sigmaz(),pulse], [sigmax(),pulse]], tlist=tlist)
print(H(10))
# initial state
psi0 = basis(2, 0)
result = mesolve(H, psi0, tlist, c_ops, [])
result2 = mesolve(H2, psi0, tlist, c_ops, [])
n = destroy(2)
a = expect(sigmaz(), result.states)
b = expect(sigmaz(), result2.states)
plt.plot(tlist, a, label = "a")
plt.plot(tlist, b, label = "b")
plt.legend()
plt.show()







