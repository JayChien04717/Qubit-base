import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
from qutip import *


#in MHz
g = 9.14  # coupling strength
g1 = 0.35   # relaxation rate
g2 = 0.035       # dephasing rate
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



def gaussian(x: np.array, peak_x: float, sigma: float):
    """
    Gaussian pulse generating function

    Parameters
    ----------
    x : np.array
        Wave event timeline, start from 0 is demanded.
    peak_x : float
        Position of Gaussian pulse peak.
    sigma : float
        Standard deviation.

    Returns
    -------
    np.array
        Output wave values.

    """
    return np.exp((-1 * (x - peak_x)**2) / (2 * sigma**2))

def gaussian_square(
        x: np.array, first_peak_x: float, flat: float, sigma: float
        ):
    """
    Unit square pulse with Gaussian edges.

    Parameters
    ----------
    x : np.array
        Wave event timeline, start from 0 is demanded.
    first_peak_x : float
        Position of the first Gaussian edge peak.
    flat : float
        Time span of flat top (1s).
    sigma : float
        Standard deviation of Gaussian edge.

    Returns
    -------
    np.array
        Output wave values.

    """
    return np.concatenate((
        gaussian(x[x <= first_peak_x], first_peak_x, sigma),
        np.ones(
            x[np.logical_and(x > first_peak_x, x <= first_peak_x + flat)].size
            ),
        gaussian(
            x[x > first_peak_x + flat], x[x > first_peak_x + flat][0], sigma
            )
        ))

h1 = 0.5*tensor(sigmaz(), qeye(2))
h2 = 0.5*tensor(qeye(2), sigmaz())
couple = g*(tensor(sigmap(), sigmam())+tensor(sigmam(), sigmap()))
H = h1+h2+couple
t_period = 1/g
t = np.linspace(0,20*t_period, 1001)
#########################################################
drive = np.sin(2*np.pi*0.5*g*t)
wave = 150*gaussian_square(t, 7*t_period,3*t_period, 0.1)*drive
# wave2 = 
z_gate = tensor(qeye(2), sigmaz())
H_t = QobjEvo([H, [z_gate, wave]],tlist=t)

psi0 = tensor(basis(2,1), basis(2,0))
result = mesolve(H_t, psi0, t, c_ops,[], options=Options(num_cpus=5))

population = sigmap()*sigmam() 
a = expect(tensor(sigmaz(), qeye(2)), result.states)
b = expect(tensor(qeye(2), sigmaz()), result.states)
plt.plot(a, label="qa")
plt.plot(b, label="qb")
plt.xlabel('us')
plt.ylabel('population')
plt.legend()
plt.show()

