// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Utilities
#include "idealreactors/utilities/Utilities"
#include "utilities/ropa/OnTheFlyROPA.h"
#include "utilities/cema/OnTheFlyCEMA.h"
#include "utilities/ontheflypostprocessing/OnTheFlyPostProcessing.h"
#include "utilities/kineticsmodifier/KineticsModifier.h"

// PolimiSoot Analyzer
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// Batch reactor
#include "idealreactors/batch/BatchReactor"

class BatchReactor
{
	public:

	BatchReactor(std::string path_kinetics);
	
	unsigned int max_number_allowed_species = 100000;
	OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML;
	OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML;

	double T, P_Pa;
	OpenSMOKE::OpenSMOKEVectorDouble omega;
	double tEnd;
    double tStart = 0.; // default 0
	double volume = 1.; // default value [1 m3]
	OpenSMOKE::BatchReactor_Type type;

	unsigned int state_variables = 0;
	bool temperature_assigned = false;
	bool pressure_assigned = false;
	bool density_assigned = false;
		
	// Temperature
	void SetTemperature(double value, std::string units);

	// Pressure
	void SetPressure(double value, std::string units);

	// Density
	double rho;
	void SetDensity(double value, std::string units);

	void CeckStatusOfGasMixture() 
	{
		if (state_variables != 2)
			OpenSMOKE::FatalErrorMessage("The status of a gas mixture requires any 2 (and only 2) among: @Temperature, @Pressure and @Density");
	}
	
	// Composition
	void SetInitialComposition(std::string initial_composition_type, std::vector<std::string> names, std::vector<double> values);

	void SetInitialComposition(std::string initial_composition_type, std::vector<double> equivalence_ratios, 
		std::string fuel_composition_type,std::vector<std::string> names_fuel, std::vector<double> values_fuel,
		std::string oxidizer_composition_type, std::vector<std::string> names_oxidizer, std::vector<double> values_oxidizer);

	// Read end time
	void SetEndTime(double value, std::string units);

	// Read start time
	// I know that this is code repetition but
	// Keep things clear
    void SetStartTime(double value, std::string units);

	// Read volume
	void SetVolume(double value, std::string units);
        
	// Read exchange area
	double exchange_area = 0.;
	void SetExchangeArea(double value, std::string units);

	// Read global thermal exchange coefficient
	double global_thermal_exchange_coefficient = 0.;
	void Set_global_thermal_exchange_coefficient(double value, std::string units);

	// Environment temperature
	double T_environment = T;
	void SetEnvironmentTemperature(double value, std::string units);

	//Type
	void SetType(std::string value);

	// Options
	OpenSMOKE::BatchReactor_Options* batch_options = new OpenSMOKE::BatchReactor_Options();

	// ODE Parameters
	OpenSMOKE::ODE_Parameters*	ode_parameters = new OpenSMOKE::ODE_Parameters();

	// Sensitivity Options
	OpenSMOKE::SensitivityAnalysis_Options* sensitivity_options;

	// On the fly ROPA
	OpenSMOKE::OnTheFlyROPA* onTheFlyROPA = new OpenSMOKE::OnTheFlyROPA(*thermodynamicsMapXML, *kineticsMapXML);

	// On the fly CEMA
	OpenSMOKE::OnTheFlyCEMA* onTheFlyCEMA = new OpenSMOKE::OnTheFlyCEMA(*thermodynamicsMapXML, *kineticsMapXML, batch_options->output_path());

	// On the fly PostProcessing
	OpenSMOKE::OnTheFlyPostProcessing* on_the_fly_post_processing = new OpenSMOKE::OnTheFlyPostProcessing(*thermodynamicsMapXML, *kineticsMapXML, batch_options->output_path());

	// Ignition Delay Times
	OpenSMOKE::IgnitionDelayTimes_Analyzer*	idt = new OpenSMOKE::IgnitionDelayTimes_Analyzer();
	
	OpenSMOKE::BatchReactor_VolumeProfile* batchreactor_volumeprofile;

	// Polimi soot
	OpenSMOKE::PolimiSoot_Analyzer* polimi_soot = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML);

	// ------------------------------------------------------------------------------------------- //
	//                              Non parametric analysis                                        //
	// ------------------------------------------------------------------------------------------- //
	std::vector<double> tempo;
	std::vector<double> temperatura;
	std::vector<double> pressione;
	std::vector<OpenSMOKE::OpenSMOKEVectorDouble> frazioni_molari;
	std::vector<OpenSMOKE::OpenSMOKEVectorDouble> frazioni_massive;
	
	void Solve();

	const std::vector<double>& GetTempo() const {return tempo;};
	const std::vector<double>& GetTemperatura() const {return temperatura;};
	const std::vector<double>& GetPressione() const {return pressione;};

	std::vector<std::vector<double>> GetMoles(std::vector<std::string> species_names)
	{
		std::vector<std::vector<double>> selected_molefractions(species_names.size(), std::vector<double>(tempo.size()));
		for(unsigned int i = 0; i < species_names.size(); i++)
		{
			unsigned int j = thermodynamicsMapXML->IndexOfSpecies(species_names[i]);
			for(unsigned int k = 0; k < tempo.size(); k++)
				selected_molefractions[i][k] = frazioni_molari[k][j];
		}
		return selected_molefractions;
	}
	
	//std::vector<double> GetMasses(std::vector<std::string> species_names)
	//{
	//	return frazioni_massive;
	//}
};

#include "OpenSMOKE_BatchReactor.hpp"
