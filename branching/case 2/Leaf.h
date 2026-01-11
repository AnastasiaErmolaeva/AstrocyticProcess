#ifndef _Leaf_H
#define _Leaf_H 

#include "Defs.h"
#include "Cell.h"
#include "Impulse.h"

class Leaf :public Cell
{
private:	//types

	struct connectedAstroParts
	{
		int nConnectedAstroParts;
		iVector partIndex;
	};

	struct connectedLeafs
	{
		int nConnectedLeafs;
		iVector leafIndex;
	};

private:	//variables

	vector<Cell*> leafArr;
	Cell * astrocyte;
	connectedAstroParts connAstroParts;
	connectedLeafs connLeafs;

	Impulse		impulse;

	FP	p_old,
		p;
	FP	q_old,
		q;
	FP	qER_old,
		qER;
	FP	z_old,
		z;
	FP  n_old,
		n;
	FP  n2_old,
		n2;
	FP  IP3_old,
		IP3;
	FP  Ca_old,
		Ca;
	FP  Na_old,
		Na;
	FP  h_old,
		h;

	FP	timeStep;
	FP  simDur;

	FP  Rcell;
	FP  Rcell2;
	FP  g_Ca;
	FP  V_m;
	FP  Ca_ext;
	FP  A_noise_leaf;

	FP V_cell;
	FP V_cell2;
	FP S_cell;
	FP S_cell2;
	FP S_VF;

	fVector k;
	fVector k_old;

	int my_ind;

	const FP rungeKuttCoeffMult[4] = { 0., 0.5, 0.5, 1. };
	const FP rungeSumKuttCoeffMult[4] = { 1., 2., 2., 1. };

public:		//functions

	void init(const int& ind, Cell* astrocyte, const vector<Cell*>& leafArr);
	void init(const vector<Cell*>& tailNet) { cout << "ERROR" << endl; system("pause"); }

	void rungeKuttStep(const int& iStep, const FP& timeAstro);
	void rungeKuttStep(const int& iPart, const int& iStep, const FP& timeAstro) { cout << "ERROR" << endl; system("pause"); }

	void saveState();
	void saveState(const int& iPart) { cout << "ERROR" << endl; system("pause"); }

	const FP& getCoeff(const int& iPart, const int& iCoeff) { return k_old[iCoeff]; }

	const FP& getP(const int& iPart) { return my_inf; }
	const FP& getQ(const int& iPart) { return my_inf; }
	const FP& getZ(const int& iPart) { return my_inf; }
	const FP& getN(const int& iPart) { return my_inf; }

	const FP& getP_old(const int& iPart) { return my_inf; }
	const FP& getQ_old(const int& iPart) { return my_inf; }
	const FP& getZ_old(const int& iPart) { return my_inf; }
	const FP& getN_old(const int& iPart) { return my_inf; }

	const FP& getCa(const int& iPart) { return Ca; }
	const FP& getNa(const int& iPart) { return Na; }
	const FP& getH(const int& iPart) { return h; }
	const FP& getN2(const int& iPart) { return n2; }

	const FP& getCa_old(const int& iPart) { return Ca_old; }
	const FP& getNa_old(const int& iPart) { return Na_old; }
	const FP& getH_old(const int& iPart) { return h_old; }
	const FP& getN2_old(const int& iPart) { return n2_old; }

private:	//functions

	void initConnectedAstroParts();

	void initConnectedLeaf();

	void constrainPrms();

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

	inline FP Noise2(const FP& nold, const FP& qold)
	{
		return -10000.*S_cell2*g_Ca*nold*(V_m - E_Ca(qold)) / (2.*96485.*V_cell2);
	}

	inline FP noise_leaf()
	{
		FP x = 0.;
		int N_cycle = 10;
		for (int i = 0; i < N_cycle; ++i)
		{
			x += FP(rand()) / FP(RAND_MAX);
		}
		x = (x - 0.5 * FP(N_cycle));
		x = x / sqrt(FP(N_cycle) / 12.);
		x = A_noise_leaf * x;
		return x;
	}

	inline FP h_inf(const FP& Caold, const FP& Naold)
	{
		return (1. - (1. / (1. + pow(Caold / KCa, HCa)))*(1. / (1. + pow(KNa / Naold, HNa))));
	}

	inline FP tau_h(const FP& Caold)
	{
		return (0.25 + tau_o / (1. + pow(Caold / Ktau, Htau)));
	}

};

#endif
