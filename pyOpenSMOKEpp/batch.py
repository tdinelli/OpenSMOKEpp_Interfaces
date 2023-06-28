import sys
sys.path.append('/Users/tdinelli/Documents/GitHub/lumpingSMOKE/source/opensmoke_interface/build')
from OpenSMOKEppInterface import BatchReactor
import numpy as np

class OS_BatchReactor:
    
	def __init__(self, kinetic_folder):
		batch = BatchReactor(kinetic_folder)
		self.batch = batch

	def InitialComposition(self, composition):
		composition_type = composition['type']
		names = composition['names']
		values = composition['values']

		self.batch.SetInitialComposition(composition_type, names, values)

	def Temperature(self, value, units):
		self.batch.SetTemperature(value, units)

	def Pressure(self, value, units):
		self.batch.SetPressure(value, units)

	def Density(self, value, units):
		self.batch.SetDensity(value, units)

	def EndTime(self, value, units):
		self.batch.SetEndTime(value, units)
	
	def StartTime(self, value, units):
		self.batch.SetStartTime(value, units)
	
	def Volume(self, value, units):
		self.batch.SetVolume(value, units)
	
	def ExchangeArea(self, value, units):
		self.batch.SetExchangeArea(value, units)

	def global_thermal_exchange_coefficient(self, value, units):
		self.batch.Set_global_thermal_exchange_coefficient(value, units)

	def EnvironmentTemperature(self, value, units):
		self.batch.SetEnvironmentTemperature(value, units)

	def Type(self, value):
		self.batch.SetType(value)

	def solve(self):
		self.batch.Solve()

	def tempo(self):
		return self.batch.tempo
	
	def temperature(self):
		return self.batch.temperature
	
	def pressure(self):
		return self.batch.pressure
	
	def mole_fractions(self, name):
		return self.batch.mole_fractions(name)
	
	def masses(self, name):
		return self.batch.mass_fractions(name)

