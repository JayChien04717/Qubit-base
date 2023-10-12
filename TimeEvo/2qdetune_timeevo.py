import numpy as np
from numpy import sqrt
import qutip as qt 
import matplotlib.pyplot as plt
from qubit import *
from op_funtion import *

import time as T

#in GHz
g = 0.75 # coupling strength
g1 = 1.15e-3    # relaxation rate
g2 = 10e-3     # dephasing rate
n_th = 10e-3      # bath temperature
c_ops = []
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



tmon1 = TunableTransmon2(13, 0.229, 0, 2, 0).hamiltonium2()
tmon2 = TunableTransmon2(12, 0.25, 0, 2, 0).hamiltonium2()

h1 = 0.5*tensor(sigmaz(), qeye(2))
h2 = 0.5*tensor(qeye(2), sigmaz())

delta = np.linspace(-10, 10, 501)
#couplint rate g = 2pi*f --> f = 1/t =  g/2pi, t =  2pi/g
time = np.linspace(0, 2*np.pi/g,101) 

z = np.vstack([np.zeros(len(time))]*len(delta))

psi0 = tensor(basis(2,1), basis(2,0))
population = sigmap()*sigmam() 



for m,n in enumerate(delta):
    H = QobjEvo([h1+h2,
    [g*(tensor(sigmap(), sigmam())),  np.exp(1j*time*n)],
    [g*(tensor(sigmam(), sigmap())), np.exp(-1j*time*n)]], 
    tlist=time)
    result = mesolve(H, psi0, time, c_ops, [], options=Options(nsteps=1e6))
    z[m,:] = expect(tensor(population, qeye(2)), result.states)

fig, ax = plt.subplots()
c = ax.imshow(z, aspect = "auto",extent = [time.min(), time.max(),delta.min(), delta.max()])
fig.colorbar(c, ax=ax)
plt.show()


