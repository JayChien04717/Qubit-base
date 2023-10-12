import numpy as np

h = 6.62607015e-34
q = 1.602e-19

def EC2C_cal():
    EC = float(input("input the EC frequency(MHz):"))
    C = q**2/(2*EC*1e6*h)
    print(f"capacitor is {C*1e15:.3f} fF") 

def C2EC_cal():
    C = float(input("input the C values(fF):"))
    EC = q**2/(2*C*1e-15*h)
    print(f"EC is {EC*1e-6:.3f} MHz") 


def EL2L_cal():
    EL = float(input("input the EL frequency(GHz):"))
    L = h/(EL*1e9*16*(np.pi*q)**2)
    print(f"inductance is {L*1e9:.2f} nH") 

def L2EL_cal():
    L = float(input("input the L values(nH):"))
    EL = h/(L*1e-9*16*(np.pi*q)**2)
    print(f"EL is {EL*1e-9:.2f} GHz") 

select = input("select the mode, 1 = EC2C, 2 = C2EC, , 3 = EL2L, 4 = L2EL: ", )

if select == str(1):
    EC2C_cal()
elif select == str(2):
    C2EC_cal()
elif select == str(3):
    EL2L_cal()
elif select == str(4):
    L2EL_cal()

