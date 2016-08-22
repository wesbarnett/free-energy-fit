
#ifndef FREEENERGYFIT_H
#define FREENERGYFIT_H

#include "LeastSquaresFit.h"
using namespace std;

class FreeEnergyFit : public LeastSquaresFit
{
    private:
        double CalcFit(vector <double> &lambda, double T);
        double ddchi2(double T, int i);
        double T0;
    public:
        // Takes initial guesses for alpha, beta, and gamma as arguments. T0 is
        // reference temperature. T is the set of temperatures corresponding
        // with the free energy data G. max_iter is the maximum number of
        // iterations to perform. stepsize is the size of step in the
        // iterations.
        FreeEnergyFit(vector <double> &lambda, double T0, vector <double> &T, vector <double> &G, int max_iter, double stepsize, double tol=1.0e-6);
        double GetTSfit(double T);
        double GetHfit(double T);
        double GetT0();
};

#endif
