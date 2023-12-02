import numpy as np
from numpy import ndarray
from qutip import *
from numpy import pi, sqrt
import matplotlib.pyplot as plt

class QubitSimulator:
    def __init__(self, number_of_qubit, couple_strength, relax=0.1, dephasing=0.1):
        '''
        In MHz unit
        '''
        self.freq = 4e9
        self.n_qubit = number_of_qubit
        self.g = couple_strength
        self.relax = relax
        self.dephas= dephasing
    def multipleQ_matrix_list(self):
        """ Function to create the multiple qubit collapse matrix.
        For 3 qubit, there has 3 list and each list conatin (2**3)x(2**3)
        matrix

        Returns
        -------
        sigmam_list: list
            The list with the sigma_minus operators in the 2^N Hilbert space. Each element in the list is 2^N*2^N dimentions
        sigmaz_list: list
            The list with the sigma_z operators in the 2^N Hilbert space. Each element in the list is 2^N*2^N dimentions
        """
        sigmam_list = [[0]]*self.n_qubit
        sigmaz_list = [[0]]*self.n_qubit

        for i in range(number_of_qubit):
            sigmam_list[i]=tensor([qeye(2)]*i+[sigmam()]+[qeye(2)]*(number_of_qubit-1-i))
            sigmaz_list[i]=tensor([qeye(2)]*i+[sigmaz()]+[qeye(2)]*(number_of_qubit-1-i))
        return sigmam_list, sigmaz_list

    def multipleQ_cops(self):
        ''' Function to give collapse operators

        Returns
        -------
        c_ops_list: list
            The list with the collapse operators in the 2^N Hilbert space.
        '''
        sigmam_list, sigmaz_list = self.multipleQ_matrix_list()
        n_th = 10e-3    # bath temperature
        c_ops=[]
        for i in range(self.n_qubit):
            c_ops.append(sqrt(self.relax * (1+n_th)) *sigmam_list[i])
            c_ops.append(sqrt(self.relax * n_th) * sigmam_list[i].dag())
            c_ops.append(sqrt(self.dephas) * sigmaz_list[i])
        return c_ops

    def couple(self):
        couple_list = [0]*(self.n_qubit-1)
        for i in range(self.n_qubit-1):
            couple_list[i] = self.g*(
                tensor([qeye(2)]*i+[sigmam()]+[sigmap()]+[qeye(2)]*(self.n_qubit-i-2))+
                tensor([qeye(2)]*i+[sigmap()]+[sigmam()]+[qeye(2)]*(self.n_qubit-i-2))
                )
        return sum(couple_list)

    def multi_qubit_hamiltonian(self):
        hamiltonian_list = [[0]]*self.n_qubit
        for i in range(self.n_qubit):
            hamiltonian_list[i] = tensor(
                [qeye(2)]*i + [0.5*self.freq*sigmaz()]+[qeye(2)]*(self.n_qubit-1-i)
            )
        return sum(hamiltonian_list)
    
    def initial_state(self, n_qubit_excited:list):
        psi0=[]
        for i in range(self.n_qubit):
            if i in n_qubit_excited:
                psi0 += [basis(2,0)]
            else:
                psi0 += [basis(2,1)]
        return tensor(psi0)

    def measure_state(self):
        population = [0]*self.n_qubit
        sz = [0]*self.n_qubit
        sx = [0]*self.n_qubit
        sy = [0]*self.n_qubit
        for i in range(self.n_qubit):
            population[i] = tensor(
                [qeye(2)]*i + [sigmap()*sigmam()]+[qeye(2)]*(self.n_qubit-1-i))
            sz[i] = tensor(
                [qeye(2)]*i + [sigmaz()]+[qeye(2)]*(self.n_qubit-1-i))
            sx[i] = tensor(
                [qeye(2)]*i + [sigmax()]+[qeye(2)]*(self.n_qubit-1-i))
            sy[i] = tensor(
                [qeye(2)]*i + [sigmax()]+[qeye(2)]*(self.n_qubit-1-i))
        return population, sz, sx, sy

    def plot(self, result, n_quit_plot:list, time:ndarray):
        pop, sz, sx, sy = self.measure_state()
        for i in n_quit_plot:
            plt.plot(time, expect(pop[i], result.states), label=f'q{i}')
        plt.xlabel("us")
        plt.legend()
        plt.show()
            


if __name__=='__main__':
    number_of_qubit = 14
    g = 0.2
    sim = QubitSimulator(number_of_qubit, g)
    H = sim.couple()+sim.multi_qubit_hamiltonian()
    c_ops = sim.multipleQ_cops()

    time_period = 1/g
    t = np.linspace(0,10*time_period, 1001)
    psi=sim.initial_state([1,3]) 
    pop, sz, sx, sy = sim.measure_state()
    # result = mesolve(H, psi, t, c_ops,[])
    result = krylovsolve(H, psi, t, krylov_dim=20, e_ops=c_ops)
    sim.plot(result, [0,1,2,3,4], t)

