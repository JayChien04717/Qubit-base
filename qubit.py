'''
This module is create the transmon base simulaor.
Author: Jay Chien
ref paper: A Quantum Engineer’s Guide to Superconducting Qubits
           Charge insensitive qubit design derived from the Cooper pair box
           Circuit quantum electrodynamic
'''

import numpy as np
from qutip import *
import op_funtion as op

class Transmon():
    def __init__(self, EJ, EC, ng) -> None:
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
        '''
        self.EJ = EJ
        self.EC = EC
        self.ng = ng
        self.ncut = ncut
        self.op = op.Operator(ncut)


    def hamiltonium(self):
        H = 4*self.EC*(self.op.n_op()-self.ng*qeye(2*self.N+1))**2 -\
            self.EJ*self.op.cos_phi_op()
        return H

    
class TunableTransmon:
    def __init__(self, EJ, EC, ng, ncut, d) -> None:
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
        self.d = d
        self.op = op.Operator(ncut)

    def hamiltonium(self,flux):
        diagonal = 4*self.EC*(self.op.n_op() - self.ng*qeye(2*self.ncut+1))**2
        off_diagonal = self.EJ*self.op.cos_phi_op()*\
            np.sqrt(np.cos(np.pi*flux)**2+self.d**2*np.sin(np.pi*flux)**2)
        H = diagonal-off_diagonal
        return H

if __name__=="__main__":
    pass