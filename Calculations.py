#NEW File
import numpy as np

L = 100

t = np.linspace(10**-3,0.02,L) #thickness from 1mm to 2 cm
w_small = np.linspace(10**-3,0.03,L) #distance between inner hole and outer radius
D = np.linspace(0.01,0.10,L) #hole diameter
W = D + w_small




A1=t*(W-D/np.sqrt(2))/2
A4=A1
A2=t*(W-D)/2
A3=A2
Aav=6/(3/A1+1/A2+1/A3+1/A4)
Abr=D*t
xaxis15 = Aav/Abr #use this for the x-axis value of graph 15
