import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
from qutip import *


#in MHz
g = 10  # coupling strength
g1 = 1.15    # relaxation rate
g2 = 0.024       # dephasing rate
n_th = 10e-3       # bath temperature

c_ops=[]
# qubit 1 collapse operators
sm1 = tensor(sigmam(), qeye(2))
sz1 = tensor(sigmaz(), qeye(2))
c_ops.append(sqrt(g1 * (1+n_th)) * sm1)
c_ops.append(sqrt(g1 * n_th) * sm1.dag())
c_ops.append(sqrt(g2) * sz1)

# qubit 2 collapse operators
sm2 = tensor(qeye(2), sigmam())
sz2 = tensor(qeye(2), sigmaz())
c_ops.append(sqrt(g1 * (1+n_th)) * sm2)
c_ops.append(sqrt(g1 * n_th) * sm2.dag())
c_ops.append(sqrt(g2) * sz2)

print(c_ops)
def const(x: np.array, lv: float):
    return lv * np.ones(x.size)

h1 = 0.5*tensor(sigmaz(), qeye(2))
h2 = 0.5*tensor(qeye(2), sigmaz())
couple = g*(tensor(sigmap(), sigmam())+tensor(sigmam(), sigmap()))
H = h1+h2+couple

t_period = 1/g
t = np.linspace(0,20*t_period, 101)
psi0 = tensor(basis(2,1), basis(2,0))
result = mesolve(couple, psi0, t, c_ops,[])

population = sigmap()*sigmam() 
a = expect(tensor(population, qeye(2)), result.states)
b = expect(tensor(qeye(2), population), result.states)
plt.plot(a, label="qa")
plt.plot(b, label="qb")
plt.xlabel('us')
plt.ylabel('population')
plt.legend()
plt.show()

