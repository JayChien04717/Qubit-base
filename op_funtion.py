import numpy as np
import scipy as sp
from numpy import ndarray
class Operator:
    def __init__(self, ncut) -> None:
        self.ncut = ncut
        self.dim = np.arange(-self.ncut, self.ncut+1,1)
    
    def n_op(self):
        """
        number present in number basis
        """
        dim = np.arange(-self.ncut, self.ncut+1,1)
        number_op = np.diag(dim)
        return number_op
    
    def exp_i_phi_op(self):
        """
        exponation(phi) operator present in number basis
        exp(i\phi)|N> = |N+1><N|
        """
        entries = np.repeat(1.0, 2*self.ncut)
        exp_op = np.diag(entries, -1) # or np.diag(-np.ones(2 * self.ncut), 1)
        return exp_op

    def cos_phi_op(self):
        """
        cos(phi) operator present in number basis
        cos(i\phi)|N> = |N+1><N|+|N-1><N| = |N+1><N|+|N><N-1|
        """
        entries = np.repeat(1.0, 2*self.ncut)
        cos_op = 0.5*(np.diag(entries, -1) + np.diag(entries, 1)) 
        # or 0.5*exp_i_phi_op()
        # cos_op += cos_op.T
        return cos_op

    def sin_phi_op(self):
        """
        sin(phi) operator present in number basis
        """
        entries = np.repeat(1.0, 2*self.ncut)
        sin_op = -1j*0.5*(np.diag(entries, -1) - np.diag(entries, 1)) 
        # or sin_op = -1j * 0.5 * self.exp_i_phi_op()
        # sin_op += sin_op.conjugate().T
        return sin_op
    
class Operator2:
    def __init__(self, dimention) -> None:
        self.dim = dimention

    def annihilation(self):
        offdiag_elements = np.sqrt(range(1, self.dim))
        return np.diagflat(offdiag_elements, 1)
    
    def creation(self):
        offdiag_elements = np.sqrt(range(1, self.dim))
        np.diagflat(offdiag_elements, 1).T
        return  self.annihilation().T #np.diagflat(offdiag_elements, 1).T

    def phi_operator(self):
        return self.creation() + self.annihilation()* self.phi_osc()/np.sqrt(2)

    def cos_phi_operator(self, alpha=1, beta=0):
        argument = alpha * self.phi_operator() + beta * np.eye(self.dim)
        native = sp.linalg.cosm(argument)
        return native

if __name__ == "__main__":
    a = Operator2(3)
    print(a.annihilation())
    print(a.creation())
    pass