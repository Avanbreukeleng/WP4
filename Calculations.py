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
tD = t/D
At = (W-D)*t
Abr = D*t

'''Mass Calculation'''
rho = 2.79E3
V   = t * (0.5 * math.pi * (W / 2) ** 2 + W * 0.8 - (math.pi * (D / 2) ** 2))
Mass = rho * V


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
kbry_006 = read("Kbry and eD for tD 0.06.txt")
kbry_08 = read("Kbry and eD for tD 0.08.txt")
kbry_010 = read("Kbry and eD for tD 0.10.txt")
kbry_012 = read("Kbry and eD for tD 0.12.txt")
kbry_015 = read("Kbry and eD for tD 0.15.txt")
kbry_060 = read("Kbry and eD for tD 0.60 and greater.txt")

A1=t*(W-D/np.sqrt(2))/2
A4=A1
A2=t*(W-D)/2
A3=A2
Aav=6/(3/A1+1/A2+1/A3+1/A4)
Abr=D*t
xaxis15 = Aav/Abr #use this for the x-axis value of graph 15

xaxis15T = xaxis15.reshape(len(xaxis15),1)
tT=t.reshape(len(t),1)
DT=D.reshape(len(D),1)
WT=W.reshape(len(W),1)
AT=At.reshape(len(At),1)
AbrT=Abr.reshape(len(Abr),1)
eDT=eD.reshape(len(eD),1)
tDT=tD.reshape(len(tD),1)
WDT=WD.reshape(len(WD),1)

values = np.hstack((tT,DT,WT,AT,AbrT,eDT,tDT,WDT,xaxis15T))
values = values[values[:, -1] < 1.36]
values = values[values[:, 5] < 4]
values = values[values[:, 5] > 0.5]
values = values[values[:, 6] > 0.055]
values = values[values[:, 7] > 1.5]
values = values[values[:, 7] < 4.9]

'''Rest of computation'''

#Stress concentration
kt = kt1_die(values[:,7])

kbr = kbry_006(values[:,5])
i = 0
for value in values[:,6]:
    if value > 0.07:
        kbr[i] = kbry_08(values[i,5])
    if value > 0.09:
        kbr[i] = kbry_010(values[i,5])
    if value > 0.11:
        kbr[i] = kbry_012(values[i,5])
    if value > 0.13:
        kbr[i] = kbry_015(values[i,5])
    if value > 0.18:
        kbr[i] = kbry_02(values[i,5])
    if value > 0.25:
        kbr[i] = kbry_03(values[i,5])
    if value > 0.35:
        kbr[i] = kbry_04(values[i,5])
    if value >=0.6:
        kbr[i] = kbry_060(values[i,5])
    i = i + 1

kty = kty(values[:,-1])
#Material properties
Fty = 365422136.54 #Pa
Ftu = 413685437.59 #Pa
density = 2800 #kg/m^3

#Allowable forces
#Py = kt * Fty * At
#Pbry = kbr * Ftu * Abr
#Pty = kty * Abr * Fty
Py = kt * Fty * values[:,3]
Pbry = kbr * Ftu * values[:,4]
Pty = kty * values[:,4] * Fty



Fa=4581.107 # Fy in Newtons
Ftr= 1471 #Fz in Newtons, taken as positive




PbryT=Pbry.reshape(len(Pbry),1)
PyT=Py.reshape(len(Py),1)
PtyT=Pty.reshape(len(Pty),1)




    
Mins = np.minimum(Pbry,Py) #finds minimum of Pbry and Pty 

Ra=Fa/Mins
Rtr=Ftr/Pty

Eq12=Ra**1.6+Rtr**1.6

MS = 1/Eq12**0.625-1

MST=MS.reshape(len(MS),1)


values=np.hstack((PyT,PbryT,PtyT,values,MST))
#values=np.hstack((PyT,PbryT,PtyT,tT,DT,WT,MST))
valuesopt=values[values[:, -1] >= 0]
valuesopt=valuesopt[valuesopt[:, -1] <= 0.5]
sortedMS=valuesopt[np.argsort(valuesopt[:, -1])]
#mass calculations, use values[:, 3 to 5]
