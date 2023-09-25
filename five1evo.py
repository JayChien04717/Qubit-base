from qutip import *
import numpy as np
from numpy import pi, sqrt
import matplotlib.pyplot as plt
#in MHz
g = 500        # coupling strength
g1 = 0.05       # relaxation rate
g2 = 0.024     # dephasing rate
n_th = 10e-3    # bath temperature

c_ops=[]
# qubit 1 collapse operators
qubitN = 16
qlist = list(np.vstack([np.zeros(1)]*qubitN))
a = [sigmam()]+ [qeye(2)]*4
aa = [sigmaz()]+ [qeye(2)]*4
b = [qeye(2)]+[sigmam()]+[qeye(2)]*3
bb = [qeye(2)]+[sigmaz()]+[qeye(2)]*3
c = [qeye(2)]*2+[sigmam()]+[qeye(2)]*2
cc = [qeye(2)]*2+[sigmaz()]+[qeye(2)]*2
d = [qeye(2)]*3+[sigmam()]+[qeye(2)]
dd = [qeye(2)]*3+[sigmaz()]+[qeye(2)]
e = [qeye(2)]*4+[sigmam()]
ee = [qeye(2)]*4+[sigmaz()]
sm1 = tensor(a)
sz1 = tensor(aa)
sm2 = tensor(b)
sz2 = tensor(bb)
sm3 = tensor(c)
sz3 = tensor(cc)
sm4 = tensor(d)
sz4 = tensor(dd)
sm5 = tensor(e)
sz5 = tensor(ee)

c_ops.append(sqrt(g1 * (1+n_th)) * sm1)
c_ops.append(sqrt(g1 * n_th) * sm1.dag())
c_ops.append(sqrt(g2) * sz1)
c_ops.append(sqrt(g1 * (1+n_th)) * sm2)
c_ops.append(sqrt(g1 * n_th) * sm2.dag())
c_ops.append(sqrt(g2) * sz2)
c_ops.append(sqrt(g1 * (1+n_th)) * sm3)
c_ops.append(sqrt(g1 * n_th) * sm3.dag())
c_ops.append(sqrt(g2) * sz3)
c_ops.append(sqrt(g1 * (1+n_th)) * sm4)
c_ops.append(sqrt(g1 * n_th) * sm4.dag())
c_ops.append(sqrt(g2) * sz4)
c_ops.append(sqrt(g1 * (1+n_th)) * sm5)
c_ops.append(sqrt(g1 * n_th) * sm5.dag())
c_ops.append(sqrt(g2) * sz5)

base = [basis(2,0)]+[basis(2,1)]*4
psi0 = tensor(base)

h1 = 0.5*tensor([sigmaz()]+ [qeye(2)]*4)
h2 = 0.5*tensor([qeye(2)]+[sigmaz()]+ [qeye(2)]*3)
h3 = 0.5*tensor([qeye(2)]*2+[sigmaz()]+ [qeye(2)]*2)
h4 = 0.5*tensor([qeye(2)]*3+[sigmaz()]+ [qeye(2)])
h5 = 0.5*tensor([qeye(2)]*4+[sigmaz()])

couple1 = g*(tensor([sigmap()]+[sigmam()]+[qeye(2)]*3) + tensor([sigmam()]+[sigmap()]+[qeye(2)]*3))
couple2 = g*(tensor([qeye(2)]+[sigmap()]+[sigmam()]+[qeye(2)]*2) + tensor([qeye(2)]+[sigmam()]+[sigmap()]+[qeye(2)]*2))
couple3 = g*(tensor([qeye(2)]*2+[sigmap()]+[sigmam()]+[qeye(2)]) + tensor([qeye(2)]*2+[sigmam()]+[sigmap()]+[qeye(2)]))
couple4 = g*(tensor([qeye(2)]*3+[sigmap()]+[sigmam()])+tensor([qeye(2)]*3+[sigmam()]+[sigmap()]))


H = h1+h2+h3+h4+h5+couple1+couple2+couple3+couple4
time_period = 1/g
t = np.linspace(0,10*time_period, 501)
result = mesolve(H, psi0, t, c_ops,[])
population = sigmap()*sigmam()
a = expect(tensor([population]+[qeye(2)]*4), result.states)
b = expect(tensor([qeye(2)]+[population]+[qeye(2)]*3), result.states)
c = expect(tensor([qeye(2)]*2+[population]+[qeye(2)]*2), result.states)
d = expect(tensor([qeye(2)]*3+[population]+[qeye(2)]), result.states)
e = expect(tensor([qeye(2)]*4+[population]), result.states)
# plt.plot(a, label = "q1")
# plt.plot(b, label = "q2")
# plt.plot(c, label = "q3")
# plt.plot(d, label = "q4")
# plt.plot(e, label = "q5")
# plt.xlabel("ns")
# plt.legend()
# plt.show()
print(qlist)