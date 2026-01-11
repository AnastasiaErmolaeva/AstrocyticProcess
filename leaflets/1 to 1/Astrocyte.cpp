#include "Astrocyte.h"

//reading parameters from "in.txt"
void Astrocyte::init(const vector<Cell*>& tailNet_in)
{
	n_Imp = 0;

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
	in >> Rcell;		
	V_cell = M_PI * Rcell * Rcell / 100.;
	S_cell = 2. * M_PI * Rcell;
	in.ignore(1000, '=');
	in >> Rcell2;		
	V_cell2 = M_PI * Rcell2 * Rcell2 / 100.;
	S_cell2 = 2. * M_PI * Rcell2 ;
	in.ignore(1000, '=');
	in >> g_Ca;
	in.ignore(1000, '=');
	in >> V_m;
	in.ignore(1000, '=');
	in >> Ca_ext;
	in.ignore(1000, '=');
	in >> A_noise;
	in.ignore(1000, '=');
	in >> A_noiseIP3;
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
	generator1 = default_random_engine(seed + 1);
	generator2 = default_random_engine(seed + 2);
	distribution = uniform_real_distribution<FP>(0., 1.);

	astroParts.resize(nAstroParts);
	initConnectedAstroParts();
	initConnectedLeaf();

	param1.resize(nAstroParts);
	param2.resize(nAstroParts);

	for (int iPart = 0; iPart < nAstroParts; ++iPart)
	{
		astroParts[iPart].c1 = 0.2367;

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
		leaf2astroPartConnections >> astroParts[iPart].connLeafs.nConnectedLeafs;
		astroParts[iPart].connLeafs.leafIndex.resize(astroParts[iPart].connLeafs.nConnectedLeafs);
		for (int iConnLeaf = 0; iConnLeaf < astroParts[iPart].connLeafs.nConnectedLeafs; ++iConnLeaf)
			leaf2astroPartConnections >> astroParts[iPart].connLeafs.leafIndex[iConnLeaf];
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

//the fourth order Runge-Kutta integration
void Astrocyte::rungeKuttStep(const int& iPart, const int& iStep, const FP& timeAstro)
{
	param1[iPart] = 1.;
	param2[iPart] = 0.;

	calcPin(iPart, iStep);
	calcQin(iPart, iStep);

	astroParts[iPart].k[0] = (astroParts[iPart].p_in + (p_1 - (astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0])) * Tr + Jplc(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1]) + param1[iPart] * noiseIP3() + param2[iPart] * noiseIP3_2()) * timeStep;
	astroParts[iPart].k[1] = (astroParts[iPart].q_in + astroParts[iPart].c1 * (Jchannel(astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0], astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1], astroParts[iPart].z_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[2], astroParts[iPart].c1) - Jpump(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1]) + Jleak(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1], astroParts[iPart].c1)) + Jin(astroParts[iPart].p_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[0]) - Jout(astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1]) + Noise(astroParts[iPart].n_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[3], astroParts[iPart].q_old + astroParts[iPart].rungeKuttCoeffMult[iStep] * astroParts[iPart].k_old[1]) ) * timeStep;
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