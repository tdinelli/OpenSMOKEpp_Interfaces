import os
import sys

if os.path.isdir("/Users/tdinelli"):
    sys.path.append(
        "/Users/tdinelli/Documents/GitHub/OpenSMOKEpp_Interfaces/build/Release/source/python")

from pyOpenSMOKE import OpenSMOKEMaps, BatchReactor

mech = OpenSMOKEMaps(
    "/Users/tdinelli/Documents/GitHub/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/13a-ignition-delay-times/kinetics",
    False)

mech.ReadMechanism()

thermo = mech.ThermodynamicsMap()
kinetics = mech.KineticsMap()

names = thermo.NamesOfSpecies()
# print(names)

reactor = BatchReactor()

reactor.SetType("NonIsothermal-ConstantVolume")
reactor.SetTemperature(700, "K")
reactor.SetPressure(20, "bar")


A_list = [kinetics.A(0), kinetics.A(0)-10000]

reactor.SetStartTime(0, "s")
reactor.SetEndTime(0.1, "s")
reactor.SetOdeOptions()
reactor.SetBatchOptions(True)


for A in A_list:
    kinetics.Set_A(0, A)
    reactor.SetKinetics(kinetics)
    reactor.SetThermodynamic(thermo)
    reactor.SetInitialComposition("MassFractions",
                                  ["N2", "O2", "NC7H16"],
                                  [7.193901e-01, 2.184260e-01, 6.218394e-02])
    reactor.SetAdditionalOptions()
    reactor.Solve()
    print(reactor.xf()[0])
