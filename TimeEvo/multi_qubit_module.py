import numpy as np
from qutip import *
from numpy import pi, sqrt
import matplotlib.pyplot as plt

class QubitSimulator:
    def __init__(self, number_of_qubit, couple_strength):
        self.n_qubit = number_of_qubit
        self.g = couple_strength
    def multipleQ_matrix_list(self):
        """ Function to create the multiple qubit collapse matrix.
        For 3 qubit, there has 3 list and each list conatin (2**3)x(2**3)
        matrix
        """
        sigmam_list = [[0]]*self.n_qubit
        sigmaz_list = [[0]]*self.n_qubit

        for i in range(number_of_qubit):
            sigmam_list[i]=tensor([qeye(2)]*i+[sigmam()]+[qeye(2)]*(number_of_qubit-1-i))
            sigmaz_list[i]=tensor([qeye(2)]*i+[sigmaz()]+[qeye(2)]*(number_of_qubit-1-i))
        return sigmam_list, sigmaz_list

    def multipleQ_cops(self):
        #in MHz
        sigmam_list, sigmaz_list = self.multipleQ_matrix_list()
        g1 = 1.15    # relaxation rate
        g2 = 0.024       # dephasing rate
        n_th = 10e-3    # bath temperature
        c_ops=[]
        for i in range(self.n_qubit):
            c_ops.append(sqrt(g1 * (1+n_th)) *sigmam_list[i])
            c_ops.append(sqrt(g1 * n_th) * sigmam_list[i].dag())
            c_ops.append(sqrt(g2) * sigmaz_list[i])
        return c_ops


    def couple(self):
        g = 1 * 2 * pi
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
                [qeye(2)]*i + [0.5*sigmaz()]+[qeye(2)]*(self.n_qubit-1-i)
            )
        return sum(hamiltonian_list)
    
    def initial_state(self, n_qubit_excited):
        psi0=[]
        for i in range(self.n_qubit):
            if i==n_qubit_excited:
                psi0 += [basis(2,0)]
            else:
                psi0 += [basis(2,1)]
        return tensor(psi0)

    def measure_state(self):
        population = [0]*self.n_qubit
        sz = 0


number_of_qubit = 2
g = 100
sim = QubitSimulator(number_of_qubit, g)
H = sim.couple()+sim.multi_qubit_hamiltonian()
c_ops = sim.multipleQ_cops()    
time_period = 1/g
t = np.linspace(0,10*time_period, 101)
base = [basis(2,0)]+[basis(2,1)]*(number_of_qubit-1)
psi0 = tensor(base)
psi=sim.initial_state(0)

result = mesolve(H, psi, t, c_ops,[])


population = sigmap()*sigmam()
a = expect(tensor([population]+[qeye(2)]*(number_of_qubit-1)), result.states)

# b = expect(tensor([qeye(2)]+[population]+[qeye(2)]*3), result.states)
# c = expect(tensor([qeye(2)]*2+[population]+[qeye(2)]*2), result.states)
# d = expect(tensor([qeye(2)]*3+[population]+[qeye(2)]), result.states
# e = expect(tensor([qeye(2)]*4+[population]), result.states)
plt.plot(a, label = "q1")
# plt.plot(b, label = "q2")
# plt.plot(c, label = "q3")
# plt.plot(d, label = "q4")
# plt.plot(e, label = "q5")
plt.xlabel("ns")
plt.legend()
plt.show()
