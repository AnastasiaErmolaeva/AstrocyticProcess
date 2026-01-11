#include "Leaf.h"
#include "Astrocyte.h"

typedef vector<Cell*> LeafArray;

void readParameters(int& nAstroParts, int& nLeafs, FP& timeStepAstro, FP& simDuration)
{
	ifstream in("input/in.txt");
	in.ignore(1000, '=');
	in >> nAstroParts;
	in.ignore(1000, '=');
	in >> nLeafs;
	string line;
	do
	{
		getline(in, line);
	} while (line.find("Simulation parameters") == string::npos);
	in.ignore(1000, '=');
	in >> timeStepAstro;
	in.ignore(1000, '=');
	in >> simDuration;

}

void initNet(Cell* astrocyte, LeafArray& leafArr, const int& nLeafs)
{
	leafArr.resize(nLeafs);
	for (int j = 0; j < nLeafs; ++j)
		leafArr[j] = new Leaf();

	for (int j = 0; j < nLeafs; ++j)
		leafArr[j]->init(j, astrocyte, leafArr);

	astrocyte->init(leafArr);
}

int main()
{
	FP simDuration, timeAstro, timeStepAstro;
	int nAstroParts, nLeafs;
	srand(int(time(NULL)));

	readParameters(nAstroParts, nLeafs, timeStepAstro, simDuration);

	Cell* astrocyte = new Astrocyte();
	LeafArray leafArr;

	initNet(astrocyte, leafArr, nLeafs);
	long long N = (long long)(ceil(long double(simDuration) / long double(timeStepAstro))), iStep;
	int i;

	vector<ofstream> of_q(nAstroParts);

	for (i = 0; i < nAstroParts; ++i)
		of_q[i].open("output/q_astroPart" + to_string(i) + ".txt");

	int nThreads = 4;
	for (iStep = 0ll; iStep < N; iStep += 1ll)
	{
		timeAstro = FP((long double)iStep * (long double)timeStepAstro);

#pragma omp parallel num_threads(nThreads)
		{
			//Runge-Kutt step 1
#pragma omp for
			for (i = 0; i < nAstroParts; ++i)
				astrocyte->rungeKuttStep(i, 0, timeAstro);
#pragma omp for
			for (i = 0; i < nLeafs; ++i)
				leafArr[i]->rungeKuttStep(0, timeAstro);

			//Runge-Kutt step 2
#pragma omp for
			for (i = 0; i < nAstroParts; ++i)
				astrocyte->rungeKuttStep(i, 1, timeAstro);
#pragma omp for
			for (i = 0; i < nLeafs; ++i)
				leafArr[i]->rungeKuttStep(1, timeAstro);

			//Runge-Kutt step 3
#pragma omp for
			for (i = 0; i < nAstroParts; ++i)
				astrocyte->rungeKuttStep(i, 2, timeAstro);
#pragma omp for
			for (i = 0; i < nLeafs; ++i)
				leafArr[i]->rungeKuttStep(2, timeAstro);

			//Runge-Kutt step 4
#pragma omp for
			for (i = 0; i < nAstroParts; ++i)
				astrocyte->rungeKuttStep(i, 3, timeAstro);
#pragma omp for
			for (i = 0; i < nLeafs; ++i)
				leafArr[i]->rungeKuttStep(3, timeAstro);

#pragma omp for
			for (i = 0; i < nAstroParts; ++i)
				astrocyte->saveState(i);
#pragma omp for
			for (i = 0; i < nLeafs; ++i)
				leafArr[i]->saveState();

#pragma omp for
			for (i = 0; i < nAstroParts; ++i)
			{
				of_q[i] << setprecision(12) << timeAstro << "\t" << astrocyte->getQ(i) << endl;
			}

#pragma omp single
			{
				if (iStep % 10000 == 0)
					cout << (N - iStep) / 10000 << endl;
			}
		}
	}

	for (i = 0; i < nAstroParts; ++i)
		of_q[i].close();

	return 0;
}