'''
This module is create the transmon base simulaor.
Author: Jay Chien
ref paper: A Quantum Engineer’s Guide to Superconducting Qubits, Charge insensitive qubit design derived from the Cooper pair box
'''

import numpy as np
import qutip as qt
import op_funtion as op

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
    def n_operator(self):
        pass
    def hamiltonium(self):
        pass

if __name__=="__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    from qutip import Qobj, about, energy_level_diagram, ket2dm, mesolve
    print(qt.sigmam())