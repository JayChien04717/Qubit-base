import numpy as np 
from numpy import pi, sqrt
import matplotlib.pyplot as plt
from pulsemodule import *
from qutip import *

def cops_term():
    g1 = 0.5    # relaxation rate
    g2 = 0.001       # dephasing rate
    n_th = 10e-3       # bath temperature
    c_ops = []
    c_ops.append(sqrt(g1 * (1+n_th)) * sigmam())
    c_ops.append(sqrt(g1 * n_th) * sigmam().dag())
    c_ops.append(sqrt(g2) * sigmaz())
    return c_ops

def pulse_obj(*args):
    pulse=[]
    for i in args:
        pulse.append(i)
    pulse = (np.array(pulse)).flatten()
    return pulse

def normalizefun(val):
    x_norm = (val-np.min(val))/(np.max(val)-np.min(val))
    return x_norm
def const(x: np.array, lv: float):
    return lv * np.ones(x.size)
def population():
    return sigmap()*sigmam()
# H = QobjEvo([0.5*sigmaz(),[0.5*sigmaz(), amplitude],[0.5*(np.cos(phase)*sigmaz()+np.sin(phase)*sigmay()), pulse]])
t = np.linspace(0,100,1001)
amplitude = np.sin(100*t)
# H = QobjEvo([0.5*sigmaz(),[0.5*sigmaz(), amplitude]],tlist=t)
H2 = QobjEvo([0.5*sigmax(),[0.5*sigmaz(),30*amplitude]], tlist=t)

t1 = np.linspace(0,20,100)
t2 = np.linspace(21, 30, 100)
tlist = np.concatenate((t1,t2))
pulse = np.concatenate((const(t1,0), const(t2,0)))
H = QobjEvo([0.5*sigmaz(), [sigmax(),pulse]], tlist=tlist)

psi0 = basis(2, 1)
result = mesolve(H, psi0, t, cops_term(), [])


a = expect(population(), result.states)
plt.plot(t, a, label = "b")
plt.legend()
plt.show()
