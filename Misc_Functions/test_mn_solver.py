import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def mn_finder(Q,Rl,Rs):
    fn = lambda Ri: Q - np.sqrt((Rl/Ri)-1) - np.sqrt((Rs/Ri)-1)
    initial_guess = min(Rs,Rl)/2
    Ri = fsolve(fn,initial_guess)
    print(Ri)
    Qr = np.sqrt((Rl/Ri)-1)
    Ql = np.sqrt((Rs/Ri)-1)
    print(Qr,Ql)
    wo = 2*np.pi*2e9
    C = (1/(wo*Ri))*(1/Q)
    print(C)
    Lload = Rl/(wo*Qr)
    print(Lload)
    Lsrc = Rs/(wo*Ql)
    print(Lsrc)  

mn_finder(3,50,35)




