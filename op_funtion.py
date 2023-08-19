import numpy as np

class Operator:
    def __init__(self, ncut) -> None:
        self.ncut = ncut
    
    def n_op(self):
        dim = np.arange(-self.ncut, self.ncut+1,1)
        chagre_basis = np.diag(dim)
        return chagre_basis
    
    def exp_i_phi_op(self):
        """
        exp(i\phi)|N> = |N+1><N|
        """
        entries = np.repeat(1.0, 2*self.ncut)
        exp_op = np.diag(entries, -1) # or np.diag(-np.ones(2 * self.ncut), 1)
        return exp_op
    def cos_phi_op(self):
        """
        cos(i\phi)|N> = |N+1><N|+|N-1><N| = |N+1><N|+|N><N-1|
        """
        entries = np.repeat(1.0, 2*self.ncut)
        cos_op = 0.5*(np.diag(entries, -1) + np.diag(entries, 1)) 
        # or 0.5*exp_i_phi_op()
        # cos_op += cos_op.T
        return cos_op
    def sin_phi_op(self):
        entries = np.repeat(1.0, 2*self.ncut)
        sin_op = -1j*0.5*(np.diag(entries, -1) - np.diag(entries, 1)) 
        # or sin_op = -1j * 0.5 * self.exp_i_phi_op()
        # sin_op += sin_op.conjugate().T
        return sin_op

if __name__ == "__main__":
    pass
    a = Operator(2)
    a = 4*(a.n_op()-0.3)**2
    print(a)
    # print(a.n_op()**2)
    # print(a.cos_phi_op())
