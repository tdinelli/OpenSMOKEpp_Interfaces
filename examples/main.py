import sys
sys.path.append('/Users/tdinelli/Documents/GitHub/lumpingSMOKE/source/opensmoke_interface/build')
from OpenSMOKEppInterface import OpenSMOKEMaps
import numpy as np
import matplotlib.pyplot as plt

maps = OpenSMOKEMaps("/Users/tdinelli/OneDrive - Politecnico di Milano/PhD/Papers/OMEs-class/mechanisms/CAI/kinetics", False)
thermo = maps.thermodynamicsMapXML
kinetics = maps.kineticsMapXML

k = []
T = np.linspace(300, 2000, 100)

for i in T:
	kinetics.SetTemperature(i)
	kinetics.KineticConstants()
	k.append(kinetics.kArrheniusModified()[1110])


x = [1000/i for i in T]
plt.plot(x, k, 'r-o')
plt.yscale('log')
plt.show()

