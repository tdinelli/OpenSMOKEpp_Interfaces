import os
import sys

if os.path.isdir("/Users/tdinelli"):
    sys.path.append(
        '/Users/tdinelli/Documents/GitHub/pyOpenSMOKEpp/build/Debug/interfaces/python')

from pyOpenSMOKE import OpenSMOKEMaps, ThermodynamicsMap_CHEMKIN
# import numpy as np
# import matplotlib.pyplot as plt

maps = OpenSMOKEMaps(
    "/Users/tdinelli/Documents/GitHub/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/01a-isothermal-constantvolume/kinetics", False)

maps.ReadMechanism()


thermo = maps.ThermodynamicsMap()
# kinetics = maps.KineticsMap()


# thermo = maps.thermodynamicsMap
# thermo = ThermodynamicsMap_CHEMKIN(maps.ThermodynamicsMap())

print("Number Of Species: ", thermo.NumberOfSpecies())
print("MWs: ", thermo.MWs())


# kinetics = maps.kineticsMap

# k = []
# T = np.linspace(300, 2000, 100)
#
# for i in T:
# 	kinetics.SetTemperature(i)
# 	kinetics.KineticConstants()
# 	k.append(kinetics.kArrheniusModified()[1110])
#
#
# x = [1000/i for i in T]
# plt.plot(x, k, 'r-o')
# plt.yscale('log')
# plt.show()
