#ifndef _ASTROCYTE_H
#define _ASTROCYTE_H 

#include "Defs.h"
#include "Cell.h"
#include "Impulse.h"

class Astrocyte :public Cell
{
private:	//types

	struct connectedAstroParts
	{
		int nConnectedAstroParts;
		iVector partIndex;
	};

	struct astroPart
	{
		connectedAstroParts connAstroParts;

		FP	p_old,
			p;
		FP	q_old,
			q;
		FP	z_old,
			z;
		FP	n_old,
			n;

		FP	p_in;
		FP	q_in;

		FP	Ca_in;
		FP	N2_in;

		FP  Rcell;
		FP  V_cell;
		FP  S_cell;
		FP  vol_ratio;

		FP	c1;		

		fVector k;
		fVector k_old;

		int connLeafInd;
		bool isLeafConnected;

		Impulse		impulse;

		const FP rungeKuttCoeffMult[4] = { 0., 0.5, 0.5, 1. };
		const FP rungeSumKuttCoeffMult[4] = { 1., 2., 2., 1. };
	};

private:	//variables

	int nAstroParts;
	vector<astroPart> astroParts;
	vector<Cell*> leafArr;

	ofstream of_IP3st;
	ifstream IP3st;

	FP	timeStep;
	FP  simDur;

	int seed;

	FP  d_IP3;
	FP  d_Ca;

	FP  g_Ca;
	FP  V_m;
	FP  Ca_ext;
	FP  A_noise;

	FP  Rcell2;
	FP  V_cell2;
	FP  S_VF;
	FP  S_cell2;

	fVector p_1;		
	int n_Imp;
	FP  impPeriod;

	fVector param;

	default_random_engine generator;
	uniform_real_distribution<FP> distribution;

public:		//functions

	Astrocyte() {};
	~Astrocyte() {};

	void init(const vector<Cell*>& leafArr);
	void init(const int& ind, Cell* astrocyte, const vector<Cell*>& leafArr) { cout << "ERROR" << endl; system("pause"); }

	void rungeKuttStep(const int& iPart, const int& iStep, const FP& timeAstro);
	void rungeKuttStep(const int& iStep, const FP& timeAstro) { cout << "ERROR" << endl; system("pause"); }

	void saveState(const int& iPart);
	void saveState() { cout << "ERROR" << endl; system("pause"); }

	const FP& getCoeff(const int& iPart, const int& iCoeff) { return astroParts[iPart].k_old[iCoeff]; }

	const FP& getP(const int& iPart) { return astroParts[iPart].p; }
	const FP& getQ(const int& iPart) { return astroParts[iPart].q; }
	const FP& getZ(const int& iPart) { return astroParts[iPart].z; }
	const FP& getN(const int& iPart) { return astroParts[iPart].n; }

	const FP& getP_old(const int& iPart) { return astroParts[iPart].p_old; }
	const FP& getQ_old(const int& iPart) { return astroParts[iPart].q_old; }
	const FP& getZ_old(const int& iPart) { return astroParts[iPart].z_old; }
	const FP& getN_old(const int& iPart) { return astroParts[iPart].n_old; }

	const FP& getCa(const int& iPart) { return my_inf; }
	const FP& getN2(const int& iPart) { return my_inf; }

	const FP& getCa_old(const int& iPart) { return my_inf; }
	const FP& getN2_old(const int& iPart) { return my_inf; }

private:	//functions

	void initConnectedAstroParts();

	void initConnectedLeaf();

	void calcPin(const int& iPart, const int& iStep);
	void calcQin(const int& iPart, const int& iStep);

	void calcCa_in(const int& iPart, const int& iStep);
	void calcN2_in(const int& iPart, const int& iStep);

	inline FP Jplc(const FP& qold)
	{
		return v4 * ((qold + (1. - alf)*k_4) / (qold + k_4));
	}

	inline FP mFun(const FP& pold)
	{
		return (pold / (pold + d1));
	}

	inline FP nFun(const FP& qold)
	{
		return (qold / (qold + d5));
	}

	inline FP Caer(const FP& qold, const FP& c1)
	{
		return (c0 - qold / (c1 * c1));
	}

	inline FP Jchannel(const FP& pold, const FP& qold, const FP& zold, const FP& c1)
	{
		return (v1 * pow(mFun(pold), 3.) * pow(nFun(qold), 3.) * pow(zold, 3.) * (Caer(qold, c1) - qold));
	}

	inline FP Jpump(const FP& qold)
	{
		return (v3 * pow(qold, 2.) / (pow(qold, 2.) + pow(k_3, 2.)));
	}

	inline FP Jleak(const FP& qold, const FP& c1)
	{
		return (v2 * (Caer(qold, c1) - qold));
	}

	inline FP Jin(const FP& pold)
	{
		return (v6 * (pow(pold, 2.) / (pow(k_2, 2.) + pow(pold, 2.))));
	}

	inline FP Jout(const FP& qold)
	{
		return k_1 * qold;
	}

	inline FP Q2(const FP& pold)
	{
		return (d2 * ((pold + d1) / (pold + d3)));
	}

	inline FP h_fun(const FP& pold, const FP& qold)
	{
		return (Q2(pold) / (Q2(pold) + qold));
	}

	inline FP Tn(const FP& pold, const FP& qold)
	{
		return (1. / (a2*(Q2(pold) + qold)));
	}

	inline FP a_m()
	{
		return 8.5 / (1. + exp(-(V_m - 8.) / 12.5));
	}

	inline FP b_m()
	{
		return 35. / (1. + exp((V_m + 74.) / 14.5));
	}

	inline FP minf()
	{
		return a_m() / (a_m() + b_m());
	}

	inline FP tauinf()
	{
		return 1. / (a_m() + b_m());
	}

	inline FP E_Ca(const FP& qold)
	{
		if (qold > 0.001)
			return 25.84*0.5*log(Ca_ext / qold);
		else
			return 127.95;
	}

	inline FP Noise(const FP& nold, const FP& qold, const FP& S_cell, const FP& V_cell)
	{
		return -10000.*S_cell*g_Ca*nold*(V_m - E_Ca(qold)) / (2.*96485.*V_cell);
	}

	inline FP Noise2(const FP& nold, const FP& qold)
	{
		return -10000.*S_cell2*g_Ca*nold*(V_m - E_Ca(qold)) / (2.*96485.*V_cell2);
	}

	inline FP noise()
	{
		FP x = 0.;
		int N_cycle = 10;
		for (int i = 0; i < N_cycle; ++i)
			x += distribution(generator);

		x = (x - 0.5 * FP(N_cycle));
		x = x / sqrt(FP(N_cycle) / 12.);
		x = A_noise * x;
		return x;
	}
};

#endif