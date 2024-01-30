import os
import sys

if os.path.isdir("/Users/tdinelli"):
    sys.path.append('/Users/tdinelli/Documents/GitHub/OpenSMOKEpp_Interfaces/build/Debug/source/python')

from pyOpenSMOKE import OpenSMOKEMaps

maps = OpenSMOKEMaps(
    "/Users/tdinelli/Documents/GitHub/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/13a-ignition-delay-times/kinetics", False)

maps.ReadMechanism()

# thermo = maps.ThermodynamicsMap()
# kinetics = maps.KineticsMap()
#
# print("Number Of Species: ", thermo.NumberOfSpecies())
# print("MWs: ", thermo.MWs())
