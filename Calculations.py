#NEW File
import numpy as np

L = 100

t = np.linspace(10**-3,0.02,L) #thickness from 1mm to 2 cm
w_small = np.linspace(10**-3,0.03,L) #distance between inner hole and outer radius
D = np.linspace(0.01,0.10,L) #hole diameter
W = D + w_small




A1=t*(W-D/np.sqrt(2))/2
