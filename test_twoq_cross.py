from qutip import *
from qubit import *
import numpy as np
import matplotlib.pyplot as plt

EJ = 13
EC = 0.24
N = 10
flux = np.linspace(-0.5, 0.5,101)
tmon = TunableTransmon(EJ, EC, 0, N,0)
tmon2 = TunableTransmon(12, EC, 0, N, 0)

z = np.zeros(len(flux))
x = np.zeros(len(flux))
# print(Qobj(op.Operator(N).n_op()).dims)
# print(tensor(qeye(tmon2.hamiltonium(0).dims[0]),tmon2.hamiltonium(0)).dims)
Nbasis=Qobj(op.Operator(N).n_op())
g=0.1
couple = g*tensor(Qobj(Nbasis), Qobj(Nbasis))

for i,j in enumerate(flux):
    h1 = tensor(tmon.hamiltonium(j), qeye(tmon.hamiltonium(j).dims[0]))
    h2 = tensor(qeye(tmon2.hamiltonium(0.2).dims[0]),tmon2.hamiltonium(0.2))
    energy = (h1+h2+couple).eigenenergies()
    z[i] = energy[1]-energy[0]
    x[i] = energy[2]-energy[0]
plt.plot(flux, z)
plt.plot(flux, x)
plt.xlim(0, 0.4)
plt.show()