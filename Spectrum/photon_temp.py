from scipy.constants import *
import numpy as np

def photon_shot(kai:float, kappa:float, dephas:float):
    """_summary_

    Args:
        kai (float): dispersive shift
        kappa (float): cavity line with
        dephas (float): qubit dephasing(in Hz unit), dephas = decoherence - 0.5*relax

    Returns:
        n_photon: _description_
    """
    n_photon = dephas*4*kappa*kai**2/(kappa**2 + 4*kai**2)
    
    return n_photon

def T_effactive(photon:float, f_res:float):
    """_summary_

    Args:
        photon (float): photon number 
        f_res (float): resonantor frequency

    Returns:
        temp: effact temperature due to photon number
    """

    temp = h*f_res/(k*np.log(1+1/photon))

    return temp