#include "Astrocyte.h"

//reading parameters from "in.txt"
void Astrocyte::init(const vector<Cell*>& tailNet_in)
{
	n_Imp = 0;
	impPeriod = 200.;

	leafArr = tailNet_in;

	ifstream in("input/in.txt");
	string line;
	in.ignore(1000, '=');
	in >> nAstroParts;
	in.ignore(1000, '=');
	in.ignore(1000, '=');
	in >> d_IP3;
	in.ignore(1000, '=');
	in >> d_Ca;
	in.ignore(1000, '=');
	in.ignore(1000, '=');
	in >> Rcell2;
	V_cell2 = M_PI * Rcell2 * Rcell2 / 100.;
	S_cell2 = 2. * M_PI * Rcell2 + 2. * M_PI * Rcell2 * Rcell2;
	in.ignore(1000, '=');
	in >> g_Ca;
	in.ignore(1000, '=');
	in >> V_m;
	in.ignore(1000, '=');
	in >> Ca_ext;
	in.ignore(1000, '=');
	in >> A_noise;
	do
	{
		getline(in, line);
	} while (line.find("Simulation parameters") == string::npos);
	in.ignore(1000, '=');
	in >> timeStep;
	in.ignore(1000, '=');
	in >> simDur;
	in.ignore(1000, '=');
	in >> seed;

	generator = default_random_engine(seed);
	distribution = uniform_real_distribution<FP>(0., 1.);

	astroParts.resize(nAstroParts);
	initConnectedAstroParts();
	initConnectedLeaf();

	param.resize(nAstroParts);

	of_IP3st.open("output/IP3st.txt");

	p_1 = fVector(nAstroParts);

	for (int iPart = 0; iPart < nAstroParts; ++iPart)
	{
		p_1[iPart] = p_1_const;

		ifstream in("input/in.txt");
		string line;
		in.ignore(1000, '=');
		in.ignore(1000, '=');
		in.ignore(1000, '=');
		in.ignore(1000, '=');
		in.ignore(1000, '=');

		in >> astroParts[iPart].Rcell;

		astroParts[iPart].V_cell = M_PI * astroParts[iPart].Rcell * astroParts[iPart].Rcell / 100.;
		astroParts[iPart].S_cell = 2. * M_PI * astroParts[iPart].Rcell;
		astroParts[iPart].vol_ratio = V_cell2 / astroParts[iPart].V_cell;

		astroParts[iPart].c1 = 0.231;

		if (iPart == 29)	
			param[iPart] = 1.;
		else 
			param[iPart] = 0.;

		astroParts[iPart].p = p0;
		astroParts[iPart].q = q0;
		astroParts[iPart].z = z0;
		astroParts[iPart].n = n0;
		astroParts[iPart].p_old = p0;
		astroParts[iPart].q_old = q0;
		astroParts[iPart].z_old = z0;
		astroParts[iPart].n_old = n0;
		astroParts[iPart].k.resize(4);
		astroParts[iPart].k_old.resize(4);
		astroParts[iPart].impulse.init(impulseAmplitude, impulseDuration, impulseStartTime, impulseFreq, impulsesQuantity, FrequencyType(freqTypeInt), AmplitudeType(ampTypeInt), impMinAmp, impMaxAmp, timeStep, iPart);
	}
}

void Astrocyte::initConnectedAstroParts()
{
	// connections between parts of one astro
	ifstream astroPartsConnections("input/astroPartsConnections.txt");
	for (int iPart = 0; iPart < nAstroParts; ++iPart)
	{
		astroPartsConnections >> astroParts[iPart].connAstroParts.nConnectedAstroParts;
		astroParts[iPart].connAstroParts.partIndex.resize(astroParts[iPart].connAstroParts.nConnectedAstroParts);
		for (int iConnAstro = 0; iConnAstro < astroParts[iPart].connAstroParts.nConnectedAstroParts; ++iConnAstro)
			astroPartsConnections >> astroParts[iPart].connAstroParts.partIndex[iConnAstro];
	}
	astroPartsConnections.close();
}

void Astrocyte::initConnectedLeaf()
{
	// connection to the leaf
	ifstream leaf2astroPartConnections("input/leaf2astroPartConnections.txt");
	int isLeafConnInt;
	for (int iPart = 0; iPart < nAstroParts; ++iPart)
	{
		leaf2astroPartConnections >> isLeafConnInt;
		astroParts[iPart].isLeafConnected = isLeafConnInt != 0;
		if (astroParts[iPart].isLeafConnected)
			leaf2astroPartConnections >> astroParts[iPart].connLeafInd;
	}
	leaf2astroPartConnections.close();
}

//IP3 diffusion
void Astrocyte::calcPin(const int& iPart, const int& iStep)
{
	astroParts[iPart].p_in = 0.;
	for (int i = 0; i < astroParts[iPart].connAstroParts.nConnectedAstroParts; ++i)
		astroParts[iPart].p_in += astroParts[astroParts[iPart].connAstroParts.partIndex[i]].p_old;

	astroParts[iPart].p_in = d_IP3 * (astroParts[iPart].p_in - astroParts[iPart].connAstroParts.nConnectedAstroParts * astroParts[iPart].p_old);
}

//Ca diffusion
void Astrocyte::calcQin(const int& iPart, const int& iStep)
{
	astroParts[iPart].q_in = 0.;

	for (int i = 0; i < astroParts[iPart].connAstroParts.nConnectedAstroParts; ++i)
		astroParts[iPart].q_in += astroParts[astroParts[iPart].connAstroParts.partIndex[i]].q_old;

	astroParts[iPart].q_in = d_Ca * (astroParts[iPart].q_in - astroParts[iPart].connAstroParts.nConnectedAstroParts * astroParts[iPart].q_old);
}

void Astrocyte::calcCa_in(const int& iPart, const int& iStep)
{
	astroParts[iPart].Ca_in = 0.;

	if (astroParts[iPart].isLeafConnected)
		astroParts[iPart].Ca_in += leafArr[astroParts[iPart].connLeafInd]->getCa_old(-1);
}

void Astrocyte::calcN2_in(const int& iPart, const int& iStep)
{
	astroParts[iPart].N2_in = 0.;

	if (astroParts[iPart].isLeafConnected)
		astroParts[iPart].N2_in += leafArr[astroParts[iPart].connLeafInd]->getN2_old(-1);
}

//the fourth order Runge-Kutta integration
void Astrocyte::rungeKuttStep(const int& iPart, const int& iStep, const FP& timeAstro)
{
	calcPin(iPart, iStep);
	calcQin(iPart, iStep);
	calcCa_in(iPart, iStep);
	calcN2_in(iPart, iStep);

	astroParts[iPart].k[0] = (astroParts[iPart].p_in + (p_1[iPart] - (astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0])) * Tr + Jplc(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1])) * timeStep;
	astroParts[iPart].k[1] = (astroParts[iPart].q_in + astroParts[iPart].c1 * (Jchannel(astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0], astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1], astroParts[iPart].z_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[2], astroParts[iPart].c1) - Jpump(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1]) + Jleak(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1], astroParts[iPart].c1)) + Jin(astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0]) - Jout(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1]) + param[iPart] * Noise(astroParts[iPart].n_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[3], astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1], astroParts[iPart].S_cell, astroParts[iPart].V_cell) + astroParts[iPart].vol_ratio * Noise2(astroParts[iPart].N2_in, astroParts[iPart].Ca_in)) * timeStep;
	astroParts[iPart].k[2] = ((h_fun(astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0], astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1]) - (astroParts[iPart].z_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[2])) / Tn(astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0], astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1])) * timeStep;
	astroParts[iPart].k[3] = ((minf() - (astroParts[iPart].n_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[3])) / tauinf() + noise()) * timeStep;

	astroParts[iPart].p += astroParts[iPart].rungeSumKuttCoeffMult[iStep] * astroParts[iPart].k[0] / 6.;
	astroParts[iPart].q += astroParts[iPart].rungeSumKuttCoeffMult[iStep] * astroParts[iPart].k[1] / 6.;
	astroParts[iPart].z += astroParts[iPart].rungeSumKuttCoeffMult[iStep] * astroParts[iPart].k[2] / 6.;
	astroParts[iPart].n += astroParts[iPart].rungeSumKuttCoeffMult[iStep] * astroParts[iPart].k[3] / 6.;

	copy(astroParts[iPart].k.begin(), astroParts[iPart].k.end(), astroParts[iPart].k_old.begin());
}

void Astrocyte::saveState(const int& iPart)
{
	astroParts[iPart].p_old = astroParts[iPart].p;
	astroParts[iPart].q_old = astroParts[iPart].q;
	astroParts[iPart].z_old = astroParts[iPart].z;
	astroParts[iPart].n_old = astroParts[iPart].n;

	astroParts[iPart].k_old = { 0.,0.,0.,0. };
}