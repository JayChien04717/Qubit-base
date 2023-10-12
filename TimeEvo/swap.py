import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
from qutip import *
g = 1.0 * 2 * pi # coupling strength
g1 = 1.15        # relaxation rate
g2 = 0.25        # dephasing rate
n_th = 0       # bath temperature

T = pi/(4*g)
H = g * (tensor(sigmax(), sigmax()) + tensor(sigmay(), sigmay()))
psi0 = tensor(basis(2,1), basis(2,0))
c_ops = []
a = destroy(4) 

# qubit 1 collapse operators
sm1 = tensor(sigmam(), qeye(2))
sz1 = tensor(sigmaz(), qeye(2))
c_ops.append(sqrt(g1 * (1+n_th)) * sm1)
c_ops.append(sqrt(g1 * n_th) * sm1.dag())
c_ops.append(sqrt(g2) * sz1)

# qubit 2 collapse operators
sm2 = tensor(qeye(2), sigmam())
sz2 = tensor(qeye(2), sigmaz())
c_ops.append(sqrt(g1 * (1+n_th)) * sm2)
c_ops.append(sqrt(g1 * n_th) * sm2.dag())
c_ops.append(sqrt(g2) * sz2)

op_basis = [[qeye(2), sigmax(), sigmay(), sigmaz()]] * 2
op_label = [["i", "x", "y", "z"]] * 2
t = np.linspace(0,1.5,101)
result = mesolve(H, psi0, t, c_ops,[])
a=destroy(2)
photon_a=tensor(sigmap()*sigmam(), qeye(2))
photon_b=tensor(qeye(2), sigmap()*sigmam())

a = expect(tensor(sigmaz(), qeye(2)), result.states)
b = expect(tensor(qeye(2), sigmaz()), result.states)
a_photon_n = expect(photon_a, result.states)
b_photon_n = expect(photon_b, result.states)
plt.plot(t,a_photon_n, label="qubit1")
plt.plot(t,b_photon_n, label="qubit2")
plt.legend()
plt.show()

