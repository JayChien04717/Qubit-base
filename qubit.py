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
    def __init__(self, EJ, EC, ng, ncut) -> None:
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
        H = 4*self.EC*(self.op.n_op()-self.ng*qeye(2*self.ncut+1))**2 -\
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

    def hamiltonium(self,flux=0):
        diagonal = 4*self.EC*(self.op.n_op() - self.ng*qeye(2*self.ncut+1))**2
        off_diagonal = self.EJ*self.op.cos_phi_op()*\
            np.sqrt(np.cos(np.pi*flux)**2+self.d**2*np.sin(np.pi*flux)**2)
        H = diagonal-off_diagonal
        return H

class Fluxnium:
    def __init__(self, EJ, EC, EL, dimention) -> None:
        self.EL = EL 
        self.EC = EC 
        self.EJ = EJ
        self.dim = dimention
        self.op = op.Operator2(dimention)
    
    def phi_osc(self) -> float:
        """
        Returns
        -------
            Returns oscillator length for the LC oscillator composed of the fluxonium
             inductance and capacitance.
        """

        return (8* self.EC / self.EL) ** 0.25  # LC oscillator length
    
    def plasma_energy(self) -> float:
        """
        Returns
        -------
            Returns the plasma oscillation frequency, sqrt(8*EL*EC).
        """
        
        return np.sqrt(8*self.EJ*self.EC) # LC plasma oscillation energy

    def hamiltonium(self, flux=0.5):
        diag_elements = [(i + 0.5) * self.plasma_energy() for i in range(self.dim)]
        lc_osc_matrix = np.diag(diag_elements)
        cos_matrix = self.op.cos_phi_operator(beta=2 * np.pi * flux, phi_osc=self.phi_osc())
        H = lc_osc_matrix - self.EJ * cos_matrix
        return H

if __name__=="__main__":
    # EJ = 9.74
    # EC = 0.99
    # EL = 0.86
    # f = Fluxnium(EJ, EC, EL, 110)

    # H = Qobj(f.hamiltonium())
    # energy = H.eigenenergies()
    # print(energy[1]-energy[0])

    def test():
        return np.random.rand(3,3)
  
    a = np.vectorize(test)
    a()
    print(np.cos(a))
    pass