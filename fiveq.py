from qutip import *
from qubit import *
import numpy as np
EJ = 10
EC = 1
EL = 1
dim = 50
couple = 0.2
flux = np.linspace(0,1,101)

# for i, flux in enumerate(flux):
#     q1 = Qobj(Fluxnium(EJ, EC, EL, dim).hamiltonium(flux))
#     q2 = Qobj(Fluxnium(EJ, EC, EL, dim).hamiltonium(flux))
#     q3 = Qobj(Fluxnium(EJ, EC, EL, dim).hamiltonium(flux))
#     q4 = Qobj(Fluxnium(EJ, EC, EL, dim).hamiltonium(flux))
#     q5 = Qobj(Fluxnium(EJ, EC, EL, dim).hamiltonium(flux))

#     h1 = tensor(q1, qeye(4*q1.dims[0]))

q1 = Qobj(Fluxnium(EJ, EC, EL, dim).hamiltonium())
h1 = tensor(q1, qeye(4*q1.dims[0]))

print(qeye(4*q1.dims[0]).shape)