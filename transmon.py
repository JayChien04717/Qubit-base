'''
This module is create the transmon base simulaor.
Author: Jay Chien
ref paper: A Quantum Engineer’s Guide to Superconducting Qubits, Charge insensitive qubit design derived from the Cooper pair box
'''

import numpy as np
import qutip as qt
import op_funtion as op
import matplotlib.pyplot as plt
class Transmon():
    def __init__(self, EJ, EC, ng, N) -> None:
        '''

        Parameter
        ----------
        EJ:
            josephson junction energy
        EC:
            capacitor energy
        ng:
            offset charge
        ncut:
            charge basis cutoff, `n = -ncut, ..., ncut`
        N:
            energy level want to calculate |N>
        '''
        self.EJ = EJ
        self.EC = EC
        self.ng = ng

        self.N = N
        self.op = op.Operator(N)


    def hamiltonium(self):
        H = 4*self.EC*(self.op.n_op() - self.ng*qt.qeye(2*self.N+1))**2 - self.EJ*self.op.cos_phi_op()
        return H

    
class TunableTransmon:
    def __init__(self, EJ, EC, ng, ncut, d, N) -> None:
        '''

        Parameter
        ----------
        EJ:
            josephson junction energy
        EC:
            capacitor energy
        ng:
            offset charge
        ncut:
            charge basis cutoff, `n = -ncut, ..., ncut`
        d:
            asymmetry between two josephton junction,
        eng_lev:
            energy level want to calculate |N>
        '''
        self.EJ = EJ
        self.EC = EC
        self.ng = ng
        self.ncut = ncut
        self.N = N
        self.d = d
        self.op = op.Operator(N)

    def hamiltonium(self,flux):
        H = 4*self.EC*(self.op.n_op() - self.ng*qt.qeye(2*self.N+1))**2 - self.EJ*self.op.cos_phi_op()*np.cos(np.pi*flux)
        return H

if __name__=="__main__":
    from qutip import *
    EJ = 13
    EC = 0.24
    N = 10
    flux = np.linspace(-0.5, 0.5,101)
    tmon = TunableTransmon(EJ, EC, 0, 0, 0, N)
    tmon2 = TunableTransmon(12, EC, 0, 0, 0, N)

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