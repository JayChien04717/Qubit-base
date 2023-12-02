import numpy as np 
from numpy import pi, sqrt
import matplotlib.pyplot as plt
from pulsemodule import *
from qutip import *
import qutip as qt
import matplotlib.colors as mcolors
"""
chalmers xmon parameter|flux pump    
f01 = 4.772            |5
Ec = 229MHz            |pump power:40MHz
EJ = 13.71GHz          | 
relaxa = 1.15MHz       |0.28MHz 
dephasing = 0.024MHz   |0.1MHz 
"""
# These terms are system collapse term
g1 = 0.5    # relaxation rate
g2 = 0.024       # dephasing rate
n_th = 10e-3       # bath temperature
c_ops = []
c_ops.append(sqrt(g1 * (1+n_th)) * sigmam()) # Jump- operator
c_ops.append(sqrt(g1 * n_th) * sigmam().dag()) # Jump+ operaotr
c_ops.append(sqrt(g2) * sigmaz()) # Pure dephasing

def const(x: np.array, lv: float):
    return lv * np.ones(x.size)



t1 = np.arange(0,15,0.01)
t2 = np.arange(15,20, 0.01)
tlist = np.concatenate((t1,t2))
pulse = np.concatenate((const(t1,2), const(t2,0)))# driving pulse
phase = 0

delta = np.linspace(-15,15,1001) #detuning

# on resonsance qubit H
H = qt.QobjEvo([0.5*sigmaz(),[0.5*(np.cos(phase)*sigmax()+np.sin(phase)*sigmay()), pulse]], tlist=tlist)

psi0 = basis(2, 1) # inistial ground state
population = sigmap()*sigmam()


z = np.vstack([np.zeros(len(tlist))]*len(delta))
for i,j in enumerate(delta):
    H = qt.QobjEvo([0.5*j*sigmaz(),[0.5*(np.cos(phase)*sigmax()+np.sin(phase)*sigmay()), pulse]], tlist=tlist)
    result = mesolve(H, psi0, tlist, c_ops, [])
    z[i,:] = expect(population, result.states)

fig, ax = plt.subplots()
c = ax.imshow(z, aspect = "auto", extent = [ pulse.min(), pulse.max(),delta.min(), delta.max()])
fig.colorbar(c, ax=ax)
plt.xlabel("us")
plt.ylabel("Detune (MHz)")
plt.show()

