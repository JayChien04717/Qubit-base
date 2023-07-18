'''
This module is create the transmon base simulaor.
Author: Jay Chien
ref paper: A Quantum Engineer’s Guide to Superconducting Qubits, Charge insensitive qubit design derived from the Cooper pair box
'''

import numpy as np
import qutip as qt


class Transmon:
    def __init__(self, EJ, EC, ng, ncut, N) -> None:
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
        self.ncut = ncut
        self.N = N

    def hamiltonium(self):
        pass






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
    from op_funtion import Operator as op
    a = op(3)
    print(op.n_op(3))