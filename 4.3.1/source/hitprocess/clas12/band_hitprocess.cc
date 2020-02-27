// G4 headers
#include "G4Poisson.hh"
#include "Randomize.hh"

// gemc headers
#include "band_hitprocess.h"

// CLHEP units
#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

// ccdb
#include <CCDB/Calibration.h>
#include <CCDB/Model/Assignment.h>
#include <CCDB/CalibrationGenerator.h>
using namespace ccdb;

static bandHitConstants initializeBANDHitConstants(int runno)
{
	// all these constants should be read from CCDB
	bandHitConstants bhc;

	// do not initialize at the beginning, only after the end of the first event,
	// with the proper run number coming from options or run table
	if(runno == -1) return bhc;

	bhc.nsector = 6;
	bhc.nlayer = 6;
	bhc.ncomp = 7;

	// database
	bhc.runNo = runno;

	bhc.date       = "2016-03-15";
	if(getenv ("CCDB_CONNECTION") != NULL)
		bhc.connection = (string) getenv("CCDB_CONNECTION");
	else
		bhc.connection = "mysql://clas12reader@clasdb.jlab.org/clas12";

	bhc.variation  = "default";
	auto_ptr<Calibration> calib(CalibrationGenerator::CreateCalibration(bhc.connection));


	vector<vector<double> > data;
	int isector, ilayer, icomp;

	cout<<"BAND:Getting effective velocities"<<endl;
	sprintf(bhc.database,"/calibration/band/effective_velocity:%d",bhc.runNo);
	data.clear(); calib->GetCalib(data,bhc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isector    = data[row][0];
		ilayer     = data[row][1];
		icomp	   = data[row][2];
		bhc.eff_vel[isector-1][ilayer-1][icomp-1] = data[row][3];
		printf("%i \t %i \t %i \t %.2f \n", isector, ilayer, icomp, bhc.eff_vel[isector-1][ilayer-1][icomp-1]);
		
	}

	cout<<"BAND:Getting TDC offsets and resolutions"<<endl;
	sprintf(bhc.database,"/calibration/band/paddle_offsets_tdc:%d",bhc.runNo);
	data.clear(); calib->GetCalib(data,bhc.database);
	for(unsigned row = 0; row < data.size(); row++)
	{
		isector    = data[row][0];
		ilayer     = data[row][1];
		icomp	   = data[row][2];
		bhc.tdc_offset[isector-1][ilayer-1][icomp-1] = data[row][3];
		bhc.tdc_resolution[isector-1][ilayer-1][icomp-1] = data[row][4];
	}

	// These are not in the CCDB
	// Fill with constant values
	cout<<"BAND:Getting MeV->ADC conversions"<<endl;
	for(isector = 0; isector < bhc.nsector; isector++) {
		for(ilayer = 0; ilayer < bhc.nlayer; ilayer++) {
			for(icomp = 0; icomp < bhc.ncomp; icomp++) {
				bhc.mev_adc[isector][ilayer][icomp] = 1000.; // channels/MeV
			}
		}
	}

	// setting voltage signal parameters
	bhc.vpar[0] = 20;  // delay, ns
	bhc.vpar[1] = 10;  // rise time, ns
	bhc.vpar[2] = 30;  // fall time, ns
	bhc.vpar[3] = 1;   // amplifier

	return bhc;
}

map<string, double> band_HitProcess :: integrateDgt(MHit* aHit, int hitn)
{
	map<string, double> dgtz;
/*
	// use Crystal ID to define IDX and IDY
	vector<identifier> identity = aHit->GetId();
	int isector    = identity[0].id;
	int ilayer     = identity[1].id;
	int icomponent = identity[2].id;

	if(aHit->isBackgroundHit == 1) {

		// background hit has all the energy in the first step. Time is also first step
		double totEdep  = aHit->GetEdep()[0];
		double stepTime = aHit->GetTime()[0];

		double charge   = totEdep*bhc.mips_charge[isector-1][ilayer-1][icomponent-1]/bhc.mips_energy[isector-1][ilayer-1][icomponent-1];

		dgtz["hitn"]      = hitn;
		dgtz["sector"]    = isector;
		dgtz["layer"]     = ilayer;
		dgtz["component"] = icomponent;
		dgtz["adc"]       = (int) (charge/bhc.fadc_LSB);
		dgtz["tdc"]       = (int) (stepTime*bhc.time_to_tdc);;

		return dgtz;
	}

	trueInfos tInfos(aHit);
	

	// initialize ADC and TDC
	int ADC = 0;
	int TDC = 8191;

	if(tInfos.eTot>0)
	{
		// adding shift and spread on time
		double time=tInfos.time+bhc.tdc_offset[isector-1][ilayer-1][icomponent-1]+G4RandGauss::shoot(0., bhc.tdc_resolution[isector-1][ilayer-1][icomponent-1]);
		TDC=int(time*bhc.time_to_tdc);
		if(TDC>bhc.tdc_max) TDC=(int)bhc.tdc_max;

		// calculate charge and amplitude
		double charge    = tInfos.eTot*bhc.mips_charge[isector-1][ilayer-1][icomponent-1]/bhc.mips_energy[isector-1][ilayer-1][icomponent-1];
		double npe_mean  = charge/bhc.gain_pc[isector-1][ilayer-1][icomponent-1];
		double npe       = G4Poisson(npe_mean);
		charge           = charge * npe/npe_mean;
		//        double amplitude = charge*bhc.gain_mv[isector-1][ilayer-1][icomponent-1]/bhc.gain_pc[isector-1][ilayer-1][icomponent-1];
		//        double fadc      = amplitude/bhc.fadc_LSB;
		ADC = (int) (charge*bhc.fadc_input_impedence/bhc.fadc_LSB/bhc.ns_per_sample);

	}

	// Status flags
	switch (bhc.status[isector-1][ilayer-1][icomponent-1])
	{
		case 0:
			break;
		case 1:
			break;
		case 3:
			ADC = TDC = 0;
			break;

		case 5:
			break;

		default:
			cout << " > Unknown FTHODO status: " << bhc.status[isector-1][ilayer-1][icomponent-1] << " for sector, layer, component "
			<< isector << ", "
			<< ilayer  << ", "
			<< icomponent << endl;
	}

	dgtz["hitn"]      = hitn;
	dgtz["sector"]    = isector;
	dgtz["layer"]     = ilayer;
	dgtz["component"] = icomponent;
	dgtz["adc"]       = ADC;
	dgtz["tdc"]       = TDC;
*/
	
	dgtz["hitn"]      = 0;
	dgtz["sector"]    = 0;
	dgtz["layer"]     = 0;
	dgtz["component"] = 0;
	dgtz["adc"]       = 0;
	dgtz["tdc"]       = 0;
	// decide if write an hit or not
	writeHit = true;

	// define conditions to reject hit
	if(rejectHitConditions) {
		writeHit = false;
	}

	return dgtz;
}

vector<identifier>  band_HitProcess :: processID(vector<identifier> id, G4Step* aStep, detector Detector)
{
	id[id.size()-1].id_sharing = 1;
	return id;
}

// - charge: returns charge/time digitized information / step
map< int, vector <double> > band_HitProcess :: chargeTime(MHit* aHit, int hitn)
{
	map< int, vector <double> >  CT;

	vector<double> hitNumbers;
	vector<double> stepIndex;
	vector<double> chargeAtElectronics;
	vector<double> timeAtElectronics;
	vector<double> identifiers;
	vector<double> hardware;
	hitNumbers.push_back(hitn);

	//////////////////////////////////////////////////////////////////
	// COMMENTED BLOCK IS COPIED FROM FT_HODO			//
	// NEED TO IMPLEMENT FOR BAND IF DIGITIZED PULSES ARE REQUIRED	//
	//////////////////////////////////////////////////////////////////

	/*

	// getting identifiers
	vector<identifier> identity = aHit->GetId();


	// use Crystal ID to define IDX and IDY
	int sector    = identity[0].id;
	int layer     = identity[1].id;
	int component = identity[2].id;
	int order = 0; // Always 0

	identifiers.push_back(sector);
	identifiers.push_back(layer);
	identifiers.push_back(component);
	identifiers.push_back(order);

	// getting hardware
	Hardware thisHardware = fthc.TT.getHardware({sector, layer, component, order});
	hardware.push_back(thisHardware.getCrate());
	hardware.push_back(thisHardware.getSlot());
	hardware.push_back(thisHardware.getChannel());

	// Adding pedestal mean and sigma into the hardware as well
	// All of these variables start from 1, therefore -1 is subtracted, e.g. sector-1
	hardware.push_back(fthc.pedestal[sector -1][layer - 1].at(component - 1));
	hardware.push_back(fthc.pedestal_rms[sector -1][layer - 1].at(component - 1));

	trueInfos tInfos(aHit);

	vector<G4ThreeVector> Lpos = aHit->GetLPos();

	vector<G4double> Edep = aHit->GetEdep();
	vector<G4double> time = aHit->GetTime();


	for (unsigned int s = 0; s < tInfos.nsteps; s++) {
		// adding shift and spread on time
		double stepTime = time[s] + fthc.time_offset[sector - 1][layer - 1][component - 1] + G4RandGauss::shoot(0., fthc.time_rms[sector - 1][layer - 1][component - 1]);

		// calculate charge and amplitude
		double stepCharge = Edep[s] * fthc.mips_charge[sector - 1][layer - 1][component - 1] / fthc.mips_energy[sector - 1][layer - 1][component - 1];
		double npe_mean = stepCharge / fthc.gain_pc[sector - 1][layer - 1][component - 1];
		double npe = G4Poisson(npe_mean);
		stepCharge = stepCharge * npe / npe_mean;
		//        double amplitude = charge*fthc.gain_mv[isector-1][ilayer-1][icomponent-1]/fthc.gain_pc[isector-1][ilayer-1][icomponent-1];
		//        double fadc      = amplitude/fthc.fadc_LSB;
		double ADC = (stepCharge * fthc.fadc_input_impedence / fthc.fadc_LSB / fthc.ns_per_sample);

		stepIndex.push_back(s);
		chargeAtElectronics.push_back(ADC);
		timeAtElectronics.push_back(stepTime);
	}

*/
	CT[0] = hitNumbers;
	CT[1] = stepIndex;
	CT[2] = chargeAtElectronics;
	CT[3] = timeAtElectronics;
	CT[4] = identifiers;
	CT[5] = hardware;

	return CT;
}

// - voltage: returns a voltage value for a given time. The inputs are:
// charge value (coming from chargeAtElectronics)
// time (coming from timeAtElectronics)
double band_HitProcess :: voltage(double charge, double time, double forTime)
{
	return PulseShape(forTime, bhc.vpar, charge, time);
}




// - electronicNoise: returns a vector of hits generated / by electronics.
vector<MHit*> band_HitProcess :: electronicNoise()
{
	vector<MHit*> noiseHits;

	// first, identify the cells that would have electronic noise
	// then instantiate hit with energy E, time T, identifier IDF:
	//
	// MHit* thisNoiseHit = new MHit(E, T, IDF, pid);

	// push to noiseHits collection:
	// noiseHits.push_back(thisNoiseHit)

	return noiseHits;
}



map< string, vector <int> >  band_HitProcess :: multiDgt(MHit* aHit, int hitn)
{
	map< string, vector <int> > MH;
	
	return MH;
}

void band_HitProcess::initWithRunNumber(int runno)
{
	if(bhc.runNo != runno) {
		cout << " > Initializing " << HCname << " digitization for run number " << runno << endl;
		bhc = initializeBANDHitConstants(runno);
		bhc.runNo = runno;
	}
}


// this static function will be loaded first thing by the executable
bandHitConstants band_HitProcess::bhc = initializeBANDHitConstants(-1);






