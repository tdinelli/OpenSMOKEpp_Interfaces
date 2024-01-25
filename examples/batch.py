import os
import sys

if os.path.isdir("/Users/tdinelli"):
    sys.path.append(
        '/Users/tdinelli/Documents/GitHub/OpenSMOKEpp_Interfaces/build/Debug/interfaces/python')

from pyOpenSMOKE import OpenSMOKEMaps, BatchReactor

mech = OpenSMOKEMaps(
    "/Users/tdinelli/Documents/GitHub/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/13a-ignition-delay-times/kinetics", True)

mech.ReadMechanism()

thermo = mech.ThermodynamicsMap()
kinetics = mech.KineticsMap()

reactor = BatchReactor()
reactor.SetType("NonIsothermal-ConstantVolume")

reactor.SetKinetics(kinetics)
reactor.SetThermodynamic(thermo)

reactor.SetTemperature(700, "K")
reactor.SetPressure(20, "bar")

reactor.SetInitialComposition("MassFractions",
                              ["N2", "O2", "NC7H16"],
                              [7.193901e-01, 2.184260e-01, 6.218394e-02])

reactor.SetStartTime(0, "s")
reactor.SetEndTime(1, "s")
reactor.SetOdeOptions()
reactor.SetBatchOptions(False)
reactor.Solve()

names = thermo.NamesOfSpecies()
xf = reactor.xf()
print(reactor.Tf())
print(reactor.Pf())
