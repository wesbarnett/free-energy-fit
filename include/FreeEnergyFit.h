
#ifndef FREEENERGYFIT_H
#define FREENERGYFIT_H

#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

extern "C" 
{
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
	void dgemm_(const char *TRANSA, const char *TRANSB, const int *M, const int *N, const int *K, double *ALPHA, double *A, const int *LDA, double *B, const int *LDB, double *BETA, double *C, const int *LDC);
}

class FreeEnergyFit
{
    private:
        vector <double> Gdata;
        vector <double> fit;
        vector <double> T;
        double calcf(double alpha, double beta, double gamma, int i);
        double ddchi2(double T, int i);
        double alpha;
        double beta;
        double gamma;
        double alpha_init;
        double beta_init;
        double gamma_init;
        double chi2;
        double T0;
        double tol;
        int max_iter;
        double stepsize;
        double CalcGfit(double alpha, double beta, double gamma, double T);
        bool converged;
    public:
        // Takes initial guesses for alpha, beta, and gamma as arguments. T0 is
        // reference temperature. T is the set of temperatures corresponding
        // with the free energy data G. max_iter is the maximum number of
        // iterations to perform. stepsize is the size of step in the
        // iterations.
        FreeEnergyFit(double alpha, double beta, double gamma, double T0, vector <double> &T, vector <double> &G, int max_iter, double stepsize, double tol=1.0e-6);
        double GetGfit(double T);
        double GetTSfit(double T);
        double GetHfit(double T);
        double GetAlpha();
        double GetBeta();
        double GetGamma();
        double GetAlphaGuess();
        double GetBetaGuess();
        double GetGammaGuess();
        double GetChi2();
        int GetMaxiter();
        double GetStepsize();
        double GetTolerance();
        bool IsConverged();
        double GetT0();
        double GetGdata(int i);
        double GetT(int i);
};

#endif
