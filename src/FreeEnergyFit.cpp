
#include "FreeEnergyFit.h"

FreeEnergyFit::FreeEnergyFit(double alpha, double beta, double gamma, double T0, vector <double> &T, vector <double> &G, int max_iter, double stepsize, double tol)
{

    this->T0 = T0;
    this->T = T;
    this->Gdata = G;
    this->chi2 = 0.0;
    this->stepsize = stepsize;
    this->max_iter = max_iter;
    this->tol = tol;
    this->converged = false;
    this->alpha_init = alpha;
    this->beta_init = beta;
    this->gamma_init = gamma;

    double f[3];
    double dx[3];
    double x[3];
    double F_prime[9];
    double d = stepsize;
    const double eps = tol;
    double sqrtf2;

    // Needed for LAPACK functions
    int N = 3;
    int N2 = 1;
    int M = 3;
    int K = 3;
    char TRANS = 'N';
    double ALPHA = 1.0;
    int LDA = 3;
    int LDB = 3;
    double BETA = 0.0;
    int LDC = 3;
    int *IPIV = new int[N];
    int LWORK = N*N;
    double *WORK = new double[LWORK];
    int INFO;

    /* Initial guesses*/
    x[0] = alpha;
    x[1] = beta;
    x[2] = gamma;

    f[0] = this->calcf(x[0], x[1], x[2], 0);
    f[1] = this->calcf(x[0], x[1], x[2], 1);
    f[2] = this->calcf(x[0], x[1], x[2], 2);

    for (int i = 0; i < max_iter; i++)
    {

        /* Forward difference approximation for derivatives */
        /* LAPACK expects the matrix as a 1-D array in column-major order, so
         * that's what we do here */
        // Col 1
        F_prime[0] = ( this->calcf(x[0] + d, x[1],     x[2],     0) - f[0] ) / d; // df0/dx0
        F_prime[1] = ( this->calcf(x[0] + d, x[1],     x[2],     1) - f[1] ) / d; // df1/dx0
        F_prime[2] = ( this->calcf(x[0] + d, x[1],     x[2],     2) - f[2] ) / d; // df2/dx0
        // Col 2
        F_prime[3] = ( this->calcf(x[0],     x[1] + d, x[2],     0) - f[0] ) / d; // df0/dx1
        F_prime[4] = ( this->calcf(x[0],     x[1] + d, x[2],     1) - f[1] ) / d; // df1/dx1
        F_prime[5] = ( this->calcf(x[0],     x[1] + d, x[2],     2) - f[2] ) / d; // df2/dx1
        // Col 3
        F_prime[6] = ( this->calcf(x[0],     x[1],     x[2] + d, 0) - f[0] ) / d; // df0/dx2
        F_prime[7] = ( this->calcf(x[0],     x[1],     x[2] + d, 1) - f[1] ) / d; // df0/dx2
        F_prime[8] = ( this->calcf(x[0],     x[1],     x[2] + d, 2) - f[2] ) / d; // df0/dx2

        dgetrf_(&N, &N, F_prime, &N, IPIV, &INFO);
        dgetri_(&N, F_prime, &N, IPIV, WORK, &LWORK, &INFO);

        /* F_prime matrix is now inverted and will now perform matrix
         * multiplication. */

        for (int j = 0; j < 9; j++)
        {
            F_prime[j] *= -1.0;
        }

        dgemm_(&TRANS, &TRANS, &M, &N2, &K, &ALPHA, F_prime, &LDA, f, &LDB, &BETA, dx, &LDC);
   
        x[0] += dx[0];
        x[1] += dx[1];
        x[2] += dx[2];

        f[0] = this->calcf(x[0], x[1], x[2], 0);
        f[1] = this->calcf(x[0], x[1], x[2], 1);
        f[2] = this->calcf(x[0], x[1], x[2], 2);

        sqrtf2 = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);

        if (sqrtf2 < eps)
        {
            this->converged = true;
            break;
        }
        
    }

    delete IPIV;
    delete WORK;
    this->alpha = x[0];
    this->beta = x[1];
    this->gamma = x[2];
    for (unsigned int i = 0; i < T.size(); i++)
    {
        chi2 += pow((this->Gdata[i] - this->CalcGfit(this->alpha, this->beta, this->gamma, this->T[i])),2);
    }
    return;

}

double FreeEnergyFit::CalcGfit(double alpha, double beta, double gamma, double T)
{

    return alpha + beta*(T-this->T0) + gamma*T*log(T/this->T0);
}

double FreeEnergyFit::GetGfit(double T)
{
    return CalcGfit(this->alpha, this->beta, this->gamma, T);
}

double FreeEnergyFit::GetTSfit(double T)
{
    return -T*( this->beta + this->gamma*( 1 + log(T/this->T0) ) );
}

double FreeEnergyFit::GetHfit(double T)
{
    return GetGfit(T) + GetTSfit(T);
}

double FreeEnergyFit::dchi2(double T, int i)
{
    // dalpha/dchi2
    if (i == 0)
    {
        return 1;
    }
    // dbeta/dchi2
    else if (i == 1)
    {
        return (T-this->T0);
    }
    // dgamma/dchi2
    else if (i == 2)
    {
        return log(T/this->T0);
    }
}

double FreeEnergyFit::calcf(double alpha, double beta, double gamma, int j)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < this->T.size(); i++)
    {
        sum += (this->Gdata[i] - this->CalcGfit(alpha, beta, gamma, this->T[i]) )*this->dchi2(T[i], j);
    }
    return sum;
}

double FreeEnergyFit::GetAlpha()
{
    return this->alpha;
}

double FreeEnergyFit::GetBeta()
{
    return this->beta;
}

double FreeEnergyFit::GetGamma()
{
    return this->gamma;
}

double FreeEnergyFit::GetAlphaGuess()
{
    return this->alpha_init;
}

double FreeEnergyFit::GetBetaGuess()
{
    return this->beta_init;
}

double FreeEnergyFit::GetGammaGuess()
{
    return this->gamma_init;
}

double FreeEnergyFit::GetChi2()
{
    return this->chi2;
}

double FreeEnergyFit::GetStepsize()
{
    return this->stepsize;
}

int FreeEnergyFit::GetMaxiter()
{
    return this->max_iter;
}

double FreeEnergyFit::GetTolerance()
{
    return this->tol;
}

double FreeEnergyFit::GetT0()
{
    return this->T0;
}

bool FreeEnergyFit::IsConverged()
{
    return this->converged;
}

double FreeEnergyFit::GetGdata(int i)
{
    return this->Gdata[i];
}

double FreeEnergyFit::GetT(int i)
{
    return this->T[i];
}

