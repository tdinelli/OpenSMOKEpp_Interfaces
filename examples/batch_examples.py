import sys
sys.path.append('/Users/tdinelli/Documents/GitHub/lumpingSMOKE/source')
from batch import OS_BatchReactor
import matplotlib.pyplot as plt
import pandas as pd

reactor = OS_BatchReactor("/Users/tdinelli/Documents/GitHub/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/01a-isothermal-constantvolume/kinetics")

reactor.Type('Isothermal-ConstantVolume')
reactor.Temperature(1000, 'K')
reactor.Pressure(101325, 'Pa')
reactor.EndTime(0.01, 's')

composition = {'type': 'MoleFractions',
               'names': ['H2', 'O2', 'N2'], 
               'values': [0.3, 0.1, 0.6]
            }

reactor.InitialComposition(composition)
reactor.solve()

t = reactor.tempo()
molefrac = reactor.mole_fractions(['OH', 'H2','N2'])

plt.plot(t, molefrac[0], '--')
plt.show()
