from qutip import *
from qubit import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

EJ = 13
EC = 0.24
N = 10
flux = np.linspace(-0.5, 0.7,101)
#|tmon, tmon2>
tmon = TunableTransmon(EJ, EC, 0, N,0) #tmon left 
tmon2 = TunableTransmon(12, EC, 0, N, 0) #tmon2 right 

z = np.zeros(len(flux))
x = np.zeros(len(flux))
y = np.zeros(len(flux))
m = np.zeros(len(flux))
# print(Qobj(op.Operator(N).n_op()).dims)
# print(tensor(qeye(tmon2.hamiltonium(0).dims[0]),tmon2.hamiltonium(0)).dims)
Nbasis=Qobj(op.Operator(N).n_op())
g=0.1
couple = g*tensor(Qobj(Nbasis), Qobj(Nbasis))
state = 6
e_state =  np.vstack([np.zeros(len(flux))]*state)
for i,j in enumerate(flux):
    h1 = tensor(tmon.hamiltonium(j), qeye(tmon.hamiltonium(j).dims[0]))
    h2 = tensor(qeye(tmon2.hamiltonium(0.2).dims[0]),tmon2.hamiltonium(0.2))
    
    energy = (h1+h2+couple).eigenenergies()
    for k in range(state):
        e_state[i,:] = energy[i]-energy[0]
        plt.plot(flux, e_state[i,:])
    # z[i] = energy[1]-energy[0]
    # x[i] = energy[2]-energy[0]
    # y[i] = energy[3]-energy[0]


# plt.plot(flux, z, label = "01")
# plt.plot(flux, x, label = "10")
# plt.plot(flux, y, label = "11")
# plt.plot(flux, m)
# plt.xlim(0, 0.7)
# plt.ylim(bottom=3.8)
# plt.legend()
plt.show()


