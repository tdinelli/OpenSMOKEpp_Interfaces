import sys
sys.path.append('/Users/tdinelli/Documents/GitHub/pyOpenSMOKEpp/pyOpenSMOKEpp')
from batch import OS_BatchReactor
import matplotlib.pyplot as plt

kinetic = "/Users/tdinelli/Documents/GitHub/OpenSMOKEppTutorials/examples/OpenSMOKEpp_BatchReactor/01a-isothermal-constantvolume/kinetics"

reactor = OS_BatchReactor(kinetic_folder=kinetic)

reactor.Type = 'Isothermal-ConstantVolume'
reactor.temperature = (1000, 'K')
reactor.pressure = (101325, 'Pa')
reactor.EndTime = (0.01, 's')

composition = {'type': 'MoleFractions',
               'species': ['H2', 'O2', 'N2'], 
               'values': [0.3, 0.1, 0.6]
            }

reactor.initial_composition = composition
reactor.solve()

t = reactor.time
molefrac = reactor.mole_fractions(['OH', 'H2','N2'])

plt.plot(t, molefrac[0], '--')
plt.show()
