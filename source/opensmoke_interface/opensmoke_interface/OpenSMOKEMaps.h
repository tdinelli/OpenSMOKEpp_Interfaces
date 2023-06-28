struct OpenSMOKEMaps
{
	OpenSMOKE::ThermodynamicsMap_CHEMKIN*	thermodynamicsMapXML;
	OpenSMOKE::KineticsMap_CHEMKIN*			kineticsMapXML;

	OpenSMOKEMaps(std::string path_kinetics, bool verbose)
	{
		boost::filesystem::path path_kinetics_ = path_kinetics;

		if (!verbose) 
			std::cout.setstate(std::ios_base::failbit);
		{
			// Read thermodynamics and kinetics maps
			boost::property_tree::ptree ptree;
    		boost::property_tree::read_xml( (path_kinetics_ / "kinetics.xml").string(), ptree );

			double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
			thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(ptree);
			kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, ptree);
			double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

			std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;
		}
		std::cout.clear();
	}
};