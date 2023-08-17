import numpy as np 
from numpy import pi, sqrt
import matplotlib.pyplot as plt
from pulsemodule import *
from qutip import *

"""
chalmers xmon parameter|flux pump    
f01 = 4.772            |5
Ec = 229MHz            |pump power:40MHz
EJ = 13.71GHz          | 
relaxa = 1.15MHz       |0.28MHz 
dephasing = 0.024MHz   |0.1MHz 
"""
def cops_term():
    g1 = 1.15    # relaxation rate
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

def const(x: np.array, lv: float):
    return lv * np.ones(x.size)

def population():
    return sigmap()*sigmam()

def h_q(f):
    return 0.5*f*sigmaz()

def q_flux(flux_envelope):
    """ 
    """
    return [0.5*sigmaz(), flux_envelope]

def q_drive(phase, pulse):
    return [0.5*(np.cos(phase)*sigmax()+np.sin(phase)*sigmay()), pulse]

t1 = np.arange(0,15,0.05)
t2 = np.arange(15,17, 0.05)
delta = np.linspace(-30,30,301)

tlist = np.concatenate((t1,t2))
flux= 20*np.sin(tlist*2*np.pi)
pulse = np.concatenate((const(t1,3), const(t2,0)))
phase1 = 0
phase2 = 0.5*np.pi


psi0 = basis(2, 1)

def plot_2D(tlist, detune):
    z = np.vstack([np.zeros(len(tlist))]*len(delta))
    for i,j in enumerate(detune):
        H2 = QobjEvo([h_q(j),q_drive(0, pulse)], tlist=tlist)
        result = mesolve(H2, psi0, tlist, cops_term(), [])
        z[i,:] = expect(population(), result.states)
    fig, ax = plt.subplots()
    ax.imshow(z,aspect = "auto")
    plt.show()
plot_2D(tlist, delta)

# z = np.vstack([np.zeros(len(tlist))]*len(delta))
# for i,j in enumerate(delta):
#     H2 = QobjEvo([h_q(j),q_drive(0, pulse)], tlist=tlist)
#     result = mesolve(H2, psi0, tlist, cops_term(), [])
#     z[i,:] = expect(population(), result.states)
# fig, ax = plt.subplots()
# ax.imshow(z,aspect = "auto")
# plt.show()
