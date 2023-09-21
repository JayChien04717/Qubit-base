from qubit import *
from qutip import *
import matplotlib.pyplot as plt

def test_TunableTransmon():
    EJ = 13
    EC = 0.29
    a = TunableTransmon(EJ, EC, 0, 31, 0)
    flux = np.linspace(-0.5,0.5,51)
    z = np.vstack([np.zeros(63)]*len(flux))
    for i, j in enumerate(flux):
        h1 = Qobj(a.hamiltonium(j))
        z[i,:] = h1.eigenenergies()
    fig, ax = plt.subplots()

    ax.plot(z[:,1]-z[:,0])
    plt.show()

def test_TunableTransmon2():
    EJ = 13
    EC = 0.29 
    b = TunableTransmon2(EJ, EC, 0, 31, 0)
    flux = np.linspace(-0.5,0.5,51)
    z2 = np.vstack([np.zeros(31)]*len(flux))
    for i, j in enumerate(flux):
        h2 = Qobj(b.hamiltonium2(j))
        z2[i,:] = h2.eigenenergies()

    fig, ax = plt.subplots()
    ax.plot(z2[:,1]-z2[:,0])
    plt.show()

def test_Fluxonium():
    #GHz
    EJ=10 #GHz
    EC=1 
    EL=1
    dim = 110
    c = Fluxnium(EJ,EC, EL, dim)
    flux = np.linspace(-0.5,0.5,101)
    z = np.vstack([np.zeros(dim)]*len(flux))
    z2 = np.vstack([np.zeros(dim)]*len(flux))
    for i, j in enumerate(flux):
        h = Qobj(c.hamiltonium2(j))
        h2 = Qobj(c.hamiltonium2(j))
        z[i,:] = h.eigenenergies()
        z2[i,:] = h2.eigenenergies()

    fig, ax = plt.subplots(1,2)
    ax[0].plot(z[:,1]-z2[:,0])
    ax[0].plot(z[:,2]-z2[:,1])
    ax[1].plot(z2[:,1]-z2[:,0])
    ax[1].plot(z2[:,2]-z2[:,1])
    
test_Fluxonium()
plt.show()

# Compare with scqubits
import scqubits as scq
fluxonium = scq.Fluxonium(
    EJ = 10,
    EC = 1,
    EL = 1,
    flux = 0.5,
    cutoff = 110
    )
fluxlist = np.linspace(-0.5,0.5,101)
energy = fluxonium.get_spectrum_vs_paramvals("flux",fluxlist , 3, False).energy_table

plt.plot(fluxlist, energy[:,1]-energy[:,0])
plt.plot(fluxlist, energy[:,2]-energy[:,1])
plt.show()