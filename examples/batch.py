import os
import sys

if os.path.isdir("/Users/tdinelli"):
    sys.path.append('/Users/tdinelli/Documents/GitHub/OpenSMOKEpp_Interfaces/build/Debug/interfaces/python')

from pyOpenSMOKE import OpenSMOKEMaps, BatchReactor

mech = OpenSMOKEMaps("/Users/tdinelli/Documents/GitHub/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/01a-isothermal-constantvolume/kinetics", False)
mech.ReadMechanism()

thermo = mech.ThermodynamicsMap()
kinetics = mech.KineticsMap()

reactor = BatchReactor()
reactor.SetType("Isothermal-ConstantPressure")

reactor.SetKinetics(kinetics)
reactor.SetThermodynamic(thermo)

reactor.SetTemperature(300, "K")
reactor.SetPressure(10, "bar")

reactor.Solve()
