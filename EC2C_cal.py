import numpy as np
def EC_cal():
    h = 6.62607015e-34
    q = 1.602e-19
    EC = float(input("inpurt the EC frequency(Hz):"))
    C = q**2/(2*EC*h)
    print(f"capacitor is {C*1e15:.2f} fF") 

def EL_cal():
    h = 6.62607015e-34
    q = 1.602e-19
    EL = float(input("inpurt the EL frequency(Hz):"))
    L = np.pi**2*h/(EL*q**2)
    print(f"capacitor is {L*1e6:.2f} nH") 

