#NEW File
import numpy as np
from scipy import interpolate
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



'''Graph Computation'''
#graph transformation
def read(file):
    x = []
    y = []
    f = open(file,"r")
    
    for line in f.readlines():
        line = line.replace(',','.')
        numbers = line.split("; ")
        numbers[1] = numbers[1].strip('\n')
        x.append(float(numbers[0] ))
        y.append(float(numbers[1]))
    f.close()
    function = interpolate.interp1d(x, y)
    return function

kty = read("Kty and Aav Abr.txt")
kt1_die = read("Kt and DW curve 1.txt")
kt2_plate = read("Kt and DW curve 2.txt")
kt3_extrusion = read("Kt and DW curve 3.txt")
kt4_bar = read("Kt and DW curve 4.txt")
kt5_hand = read("Kt and DW curve 5.txt")
kt6_alloy = read("Kt and DW curve 6.txt")
kt7_mag = read("Kt and DW curve 7.txt")
kbry_02 = read("Kbry and eD for tD 0.2.txt")
kbry_03 = read("Kbry and eD for tD 0.3.txt")
kbry_04 = read("Kbry and eD for tD 0.4.txt")
kbry_06 = read("Kbry and eD for tD 0.06.txt")
kbry_08 = read("Kbry and eD for tD 0.08.txt")
kbry_010 = read("Kbry and eD for tD 0.10.txt")
kbry_012 = read("Kbry and eD for tD 0.12.txt")
kbry_015 = read("Kbry and eD for tD 0.15.txt")


A1=t*(W-D/np.sqrt(2))/2
A4=A1
A2=t*(W-D)/2
A3=A2
Aav=6/(3/A1+1/A2+1/A3+1/A4)
Abr=D*t
xaxis15 = Aav/Abr #use this for the x-axis value of graph 15



'''Rest of computation'''

#Stress concentration
kt = 1
kbr = 1
kty = kty(xaxis15)
#Material properties
Fty = 414 * 10**6 #Pa
Ftu = 483 * 10**6 #Pa
density = 2800 #kg/m^3

#Allowable forces
Py = kt * Fty * At
Pbry = kbr * Ftu * Abr
Pty = kty * Abr * Fty



Fa=4581.107 # Fy in Newtons
Ftr= 1471 #Fz in Newtons, taken as positive




PbryT=Pbry.reshape(len(Pbry),1)
PyT=Py.reshape(len(Py),1)
PtyT=Pty.reshape(len(Pty),1)
tT=t.reshape(len(t),1)
DT=D.reshape(len(D),1)
WT=W.reshape(len(W),1)



    
Mins = np.minimum(Pbry,Py) #finds minimum of Pbry and Pty 

Ra=Fa/Mins
Rtr=Ftr/Pty

Eq12=Ra**1.6+Rtr**1.6

MS = 1/Eq12**0.625-1

MST=MS.reshape(len(MS),1)


values=np.hstack((PyT,PbryT,PtyT,tT,DT,WT,MST))
valuesopt=values[values[:, -1] >= 0]
#mass calculations
