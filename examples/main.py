import sys
sys.path.append('/Users/tdinelli/Documents/GitHub/lumpingSMOKE/source/opensmoke_interface/build')
from OpenSMOKEppInterface import KineticMap
import numpy as np
import matplotlib.pyplot as plt

km = KineticMap("/Users/tdinelli/OneDrive - Politecnico di Milano/PhD/Papers/OMEs-class/mechanisms/CAI/kinetics", False)
thermo = km.thermodynamicsMapXML
kinetics = km.kineticsMapXML

k = []
T = np.linspace(300, 2000, 100)

for i in T:
	kinetics.SetTemperature(i)
	kinetics.KineticConstants()
	k.append(kinetics.kArrheniusModified()[1110])


plt.plot(T, k, 'r-o')
plt.show()

