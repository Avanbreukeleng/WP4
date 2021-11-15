#NEW File
import numpy as np

L = 100

t_base = np.linspace(10**-3,0.02,L) #thickness from 1mm to 2 cm
w_small_base = np.linspace(10**-3,0.03,L) #distance between inner hole and outer radius
D_base = np.linspace(0.01,0.10,L) #hole diameter

t = []
D = []
w_small = []

for ti in t_base:
    for w_smalli in w_small_base:
        for Di in D_base:
            t.append(ti)
            w_small.append(w_smalli)
            D.append(Di)
t = np.array(t)
w_small = np.array(w_small)
D = np.array(D)

W = D + w_small
e = w_small + D/2
eD = e/D
WD = W/D
At = (W-D)*t
Abr = D*t
#Stress concentration
kt = 1
kbr = 1
kty = 1
#Material properties
Fty = 414 * 10**6 #Pa
Ftu = 483 * 10**6 #Pa
density = 2800 #kg/m^3

#Allowable forces
Py = kt * Fty * At
Pbry = kbr * Ftu * Abr
Pty = kty * Abr * Fty


A1=t*(W-D/np.sqrt(2))/2
A4=A1
A2=t*(W-D)/2
A3=A2
Aav=6/(3/A1+1/A2+1/A3+1/A4)
Abr=D*t
xaxis15 = Aav/Abr #use this for the x-axis value of graph 15
