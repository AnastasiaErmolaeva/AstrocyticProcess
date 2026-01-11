#include "Leaf.h"

//reading parameters from "in.txt"
void Leaf::init(const int& ind_in, Cell* astrocyte_in, const vector<Cell*>& leafArr_in)
{
	astrocyte = astrocyte_in;
	leafArr = leafArr_in;
	my_ind = ind_in;

	ifstream in("input/in.txt");
	string line;
	in.ignore(1000, '=');
	in.ignore(1000, '=');
	in.ignore(1000, '=');
	in.ignore(1000, '=');
	in.ignore(1000, '=');
	in >> Rcell;
	V_cell = M_PI * Rcell * Rcell / 100.;
	S_cell = 2. * M_PI * Rcell;
	in.ignore(1000, '=');
	in >> Rcell2;
	V_cell2 = M_PI * Rcell2 * Rcell2 / 100.;
	S_cell2 = 2. * M_PI * Rcell2 + 2. * M_PI * Rcell2*Rcell2;
	in.ignore(1000, '=');
	in >> g_Ca;
	in.ignore(1000, '=');
	in >> V_m;
	in.ignore(1000, '=');
	in >> Ca_ext;
	in.ignore(1000, '=');
	in.ignore(1000, '=');
	in >> A_noise_leaf;
	do
	{
		getline(in, line);
	} while (line.find("Simulation parameters") == string::npos);
	in.ignore(1000, '=');
	in >> timeStep;
	in.ignore(1000, '=');
	in >> simDur;

	k.resize(4);
	k_old.resize(4);

	Ca = Ca_0;
	Na = Na_0;
	h = h_0;
	n2 = n2_0;

	Ca_old = Ca_0;
	Na_old = Na_0;
	h_old = h_0;
	n2_old = n2_0;

	initConnectedAstroParts();
	initConnectedLeaf();

	if (connAstroParts.nConnectedAstroParts != 0)
		impulse.init(impulseAmplitude, impulseDuration, impulseStartTime, impulseFreq, impulsesQuantity, FrequencyType(freqTypeInt), AmplitudeType(ampTypeInt), impMinAmp, impMaxAmp, timeStep, connAstroParts.partIndex[0]);
	else
		impulse.init(impulseAmplitude, impulseDuration, impulseStartTime, impulseFreq, impulsesQuantity, FREQ_NEVER, AmplitudeType(ampTypeInt), impMinAmp, impMaxAmp, timeStep, -1);
}

void Leaf::initConnectedLeaf()
{
	// connection to the leaf
	ifstream leaf2LeafConnections("input/leafsConnections.txt");
	string line;
	for (int i = 0; i < my_ind; ++i)
		getline(leaf2LeafConnections, line);

	leaf2LeafConnections >> connLeafs.nConnectedLeafs;

	connLeafs.leafIndex.resize(connLeafs.nConnectedLeafs);
	for (int i = 0; i < connLeafs.nConnectedLeafs; ++i)
	{
		leaf2LeafConnections >> connLeafs.leafIndex[i];
	}

	leaf2LeafConnections.close();
}

void Leaf::initConnectedAstroParts()
{
	// connections to the parts of astrocyte
	ifstream leaf2astroPartConnections("input/leaf2astroPartConnections.txt");

	string line;
	int iAstroPart = -1;
	do
	{
		++iAstroPart;
		getline(leaf2astroPartConnections, line);
	} while (line.find("1\t" + to_string(my_ind)) == string::npos && iAstroPart < 1000);

	if (line.find("1\t" + to_string(my_ind)) == string::npos)
	{
		connAstroParts.nConnectedAstroParts = 0;
	}
	else
	{
		connAstroParts.nConnectedAstroParts = 1;
		connAstroParts.partIndex.resize(connAstroParts.nConnectedAstroParts);
		connAstroParts.partIndex[0] = iAstroPart;
	}
	leaf2astroPartConnections.close();
}

//the fourth order Runge-Kutta integration
void Leaf::rungeKuttStep(const int& iStep, const FP& timeAstro)
{
	k[0] = ((Ca_rest - (Ca_old + rungeKuttCoeffMult[iStep] * k_old[0])) / tau_Ca + Noise2(n2_old + rungeKuttCoeffMult[iStep] * k_old[3], Ca_old + rungeKuttCoeffMult[iStep] * k_old[0]))*timeStep;
	k[1] = ((Na_rest - (Na_old + rungeKuttCoeffMult[iStep] * k_old[1])) / tau_Na)*timeStep;
	k[2] = ((h_inf(Ca_old + rungeKuttCoeffMult[iStep] * k_old[0], Na_old + rungeKuttCoeffMult[iStep] * k_old[1]) - (h_old + rungeKuttCoeffMult[iStep] * k_old[2])) / tau_h(Ca_old + rungeKuttCoeffMult[iStep] * k_old[0]))*timeStep;
	k[3] = ((minf() - (n2_old + rungeKuttCoeffMult[iStep] * k_old[3])) / tauinf() + noise_leaf())*timeStep;

	Ca += rungeSumKuttCoeffMult[iStep] * k[0] / 6.;
	Na += rungeSumKuttCoeffMult[iStep] * k[1] / 6.;
	h += rungeSumKuttCoeffMult[iStep] * k[2] / 6.;
	n2 += rungeSumKuttCoeffMult[iStep] * k[3] / 6.;

	copy(k.begin(), k.end(), k_old.begin());

	constrainPrms();
}

void Leaf::saveState()
{
	Ca_old = Ca;
	Na_old = Na;
	h_old = h;
	n2_old = n2;

	k_old = { 0.,0.,0.,0. };
}

void Leaf::constrainPrms()
{
	if (Ca < 0.01)
	{
		Ca = 0.01;
	}
}
