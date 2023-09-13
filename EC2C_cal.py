import numpy as np

h = 6.62607015e-34
q = 1.602e-19

def EC_cal():
    EC = float(input("inpurt the EC frequency(Hz):"))
    C = q**2/(2*EC*h)
    print(f"capacitor is {C*1e15:.2f} fF") 

def EL_cal():
    EL = float(input("inpurt the EL frequency(Hz):"))
    L = np.pi**2*h/(EL*q**2)
    print(f"inductance is {L*1e6:.2f} nH") 

def capacitance_cal(target_freq):
    z_im = float(input("inpurt the imaginary part of impencdence:"))
    capacitor = 1/(2*np.pi*z_im*target_freq)
    print(f"capacitor is {capacitor*1e15:.2f} fF")
    