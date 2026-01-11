#ifndef _CELL_H
#define _CELL_H

#include "Defs.h"


class Cell
{
public:

	const FP	my_inf = INFINITY;
	const int	my_int_max = INT_MAX;

	virtual void init(const vector<Cell*>& tailNet) = 0;
	virtual void init(const int& ind, Cell* astrocyte, const vector<Cell*>& leafArr) = 0;

	virtual void rungeKuttStep(const int& iStep, const FP& timeAstro) = 0;
	virtual void rungeKuttStep(const int& iPart, const int& iStep, const FP& timeAstro) = 0;

	virtual void saveState() = 0;
	virtual void saveState(const int& iPart) = 0;

	virtual const FP& getP(const int& iPart) = 0;
	virtual const FP& getQ(const int& iPart) = 0;
	virtual const FP& getZ(const int& iPart) = 0;
	virtual const FP& getN(const int& iPart) = 0;

	virtual const FP& getP_old(const int& iPart) = 0;
	virtual const FP& getQ_old(const int& iPart) = 0;
	virtual const FP& getZ_old(const int& iPart) = 0;
	virtual const FP& getN_old(const int& iPart) = 0;

	virtual const FP& getCa(const int& iPart) = 0;
	virtual const FP& getN2(const int& iPart) = 0;

	virtual const FP& getCa_old(const int& iPart) = 0;
	virtual const FP& getN2_old(const int& iPart) = 0;

	virtual const FP& getCoeff(const int& iPart, const int& iCoeff) = 0;
};

#endif