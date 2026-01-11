#ifndef _DEFS_H
#define _DEFS_H
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>

using namespace std;

typedef double FP;
typedef vector<int> iVector;
typedef vector<FP> fVector;
typedef vector<iVector> iMatrix;
typedef vector<fVector> fMatrix;

//initial values of variables
const FP p0 = 0.16, q0 = 0.07, z0 = 0.67, n0 = 0., Ca_0 = 0.08, Na_0 = 16., h_0 = 0.1, n2_0 = 0.,

//astrocyte parameters
a2 = 0.14, d1 = 0.13, d2 = 1.049, d3 = 0.9434, d5 = 0.082,
v_b = 0.5, K_r = 1.3, K_p = 10., K_n = 0.6,
k_1 = 0.5, k_2 = 1., k_3 = 0.1, k_4 = 1.1,
c0 = 20.,
v1 = 6., v2 = 0.11, v3 = 20., v4 = 0.4, v6 = 0.2,
alf = 0.8, Tr = 0.14,
p_1 = 0.28,

//leaflet parameters
Ca_rest = 0.075, Na_rest = 16., tau_Ca = 2., tau_Na = 6.,
KCa = 0.8, HCa = 2., KNa = 5., HNa = 2.,
tau_o = 10., Ktau = 1., Htau = 1.,

//input signal parameters
impDuration = 0.1, impAmpl = 100., impStartTime = 50., impPeriod = 50.;

enum FrequencyType
{
	FREQ_POISSON = 0,
	FREQ_CONST = 1,
	FREQ_NEVER = 2
};

enum AmplitudeType
{
	AMP_CONST = 0,
	AMP_RAND_UNIFORM = 1
};

#endif