from OpenSMOKEppInterface import BatchReactor

class OS_BatchReactor:
    
	def __init__(self, kinetic_folder: str):
		batch = BatchReactor(kinetic_folder)
		self._batch = batch
		self._initial_composition = None
		self._temperature = None
		self._pressure = None
		self._density = None
		self._end_time = None
		self._start_time = None
		self._volume = None
		self._exchange_area = None
		self._global_thermal_exchange_coefficient = None
		self._environment_temperature = None
		self._reactor_type = None

	@property
	def initial_composition(self):
		return self._initial_composition

	@initial_composition.setter
	def initial_composition(self, composition):
		composition_type = composition['type']
		names = composition['species']
		values = composition['values']
		self._batch.SetInitialComposition(composition_type, names, values)
		self._initial_composition = composition

	@property
	def temperature(self):
		return self._temperature

	@temperature.setter
	def temperature(self, value):
		# value[0] is the value
		# value[1] is the UoM
		self._temperature = value[0]
		self._batch.SetTemperature(value[0], value[1])

	@property
	def pressure(self):
		return self._pressure
	
	@pressure.setter
	def pressure(self, value):
		self._pressure = value[0]
		self._batch.SetPressure(value[0], value[1])

	@property
	def density(self):
		return self._density
	
	@density.setter
	def density(self, value):
		self._density = value[0]
		self._batch.SetDensity(value[0], value[1])

	@property
	def EndTime(self):
		return self._end_time

	@EndTime.setter
	def EndTime(self, value):
		self._end_time = value[0]
		self._batch.SetEndTime(value[0], value[1])

	@property
	def StartTime(self):
		return self._start_time

	@StartTime.setter
	def StartTime(self, value):
		self._start_time = value[0]
		self._batch.SetStartTime(value[0], value[1])

	@property
	def volume(self):
		return self._volume
		
	@volume.setter
	def volume(self, value):
		self._volume = value[0]
		self._batch.SetVolume(value[0], value[1])
	
	@property
	def ExchangeArea(self):
		return self._exchange_area

	@ExchangeArea.setter
	def ExchangeArea(self, value):
		self._exchange_area = value[0]
		self._batch.SetExchangeArea(value[0], value[1])
	
	@property
	def global_thermal_exchange_coefficient(self):
		return self._global_thermal_exchange_coefficient
	
	@global_thermal_exchange_coefficient.setter
	def global_thermal_exchange_coefficient(self, value):
		self._global_thermal_exchange_coefficient = value[0]
		self._batch.Set_global_thermal_exchange_coefficient(value[0], value[1])

	@property
	def environment_temperature(self):
		return self._environment_temperature
	
	@environment_temperature.setter
	def EnvironmentTemperature(self, value):
		self._environment_temperature = value[0]
		self._batch.SetEnvironmentTemperature(value[0], value[1])

	@property
	def Type(self):
		return self._reactor_type
	
	@Type.setter
	def Type(self, value):
		self._reactor_type = value
		self._batch.SetType(value)

	def solve(self):
		self._batch.Solve()

	@property
	def time(self):
		return self._batch.time
	
	@property
	def final_temperature(self):
		return self._batch.temperature
	
	@property
	def final_pressure(self):
		return self._batch.pressure
	
	def mole_fractions(self, name):
		return self._batch.mole_fractions(name)
	
	def masses(self, name):
		return self._batch.mass_fractions(name)

