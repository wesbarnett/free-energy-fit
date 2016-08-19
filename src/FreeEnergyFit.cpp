
#include "FreeEnergyFit.h"

FreeEnergyFit::FreeEnergyFit(vector <double> &lambda, double T0, vector <double> &T, vector <double> &G, int max_iter, double stepsize, double tol)
{

    this->T0 = T0;
    this->T = T;
    this->Gdata = G;
    this->chi2 = 0.0;
    this->stepsize = stepsize;
    this->max_iter = max_iter;
    this->tol = tol;
    this->converged = false;
    for (int i = 0; i < lambda.size(); i++)
    {
        this->lambda_init.push_back(lambda[i]);
    }
    this->lambda.resize(lambda.size());

    // Needed for LAPACK functions
    int N = lambda.size();
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

    double f[N];
    double dx[N];
    vector <double> x(N);
    double F_prime[N*N];
    double d = stepsize;
    const double eps = tol;
    double sqrtf2;

    /* Initial guesses*/
    for (int i = 0; i < x.size(); i++)
    {
        x[i] = lambda[i];
    }

    for (int i = 0; i < x.size(); i++)
    {
        f[i] = this->calcf(x, i);
    }

    for (int i = 0; i < max_iter; i++)
    {

        /* Forward difference approximation for derivatives */
        /* LAPACK expects the matrix as a 1-D array in column-major order, so
         * that's what we do here */
        // Col 1
        lambda[0] = x[0] + d;
        lambda[1] = x[1];
        lambda[2] = x[2];
        F_prime[0] = ( this->calcf(lambda, 0) - f[0] ) / d; // df0/dx0
        F_prime[1] = ( this->calcf(lambda, 1) - f[1] ) / d; // df1/dx0
        F_prime[2] = ( this->calcf(lambda, 2) - f[2] ) / d; // df2/dx0
        // Col 2
        lambda[0] = x[0];
        lambda[1] = x[1] + d;
        lambda[2] = x[2];
        F_prime[3] = ( this->calcf(lambda, 0) - f[0] ) / d; // df0/dx1
        F_prime[4] = ( this->calcf(lambda, 1) - f[1] ) / d; // df1/dx1
        F_prime[5] = ( this->calcf(lambda, 2) - f[2] ) / d; // df2/dx1
        // Col 3
        lambda[0] = x[0];
        lambda[1] = x[1];
        lambda[2] = x[2] + d;
        F_prime[6] = ( this->calcf(lambda, 0) - f[0] ) / d; // df0/dx2
        F_prime[7] = ( this->calcf(lambda, 1) - f[1] ) / d; // df0/dx2
        F_prime[8] = ( this->calcf(lambda, 2) - f[2] ) / d; // df0/dx2

        dgetrf_(&N, &N, F_prime, &N, IPIV, &INFO);
        dgetri_(&N, F_prime, &N, IPIV, WORK, &LWORK, &INFO);

        /* F_prime matrix is now inverted and will now perform matrix
         * multiplication. */

        for (int j = 0; j < 9; j++)
        {
            F_prime[j] *= -1.0;
        }

        dgemm_(&TRANS, &TRANS, &M, &N2, &K, &ALPHA, F_prime, &LDA, f, &LDB, &BETA, dx, &LDC);
   
        for (int j = 0; j < x.size(); j++)
        {
            x.at(j) += dx[j];
        }

        for (int j = 0; j < x.size(); j++)
        {
            f[j] = this->calcf(x, j);
        }

        sqrtf2 = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);

        if (sqrtf2 < eps)
        {
            this->converged = true;
            break;
        }
        
    }

    delete IPIV;
    delete WORK;
    for (int i = 0; i < N; i++)
    {
        this->lambda.at(i) = x.at(i);
    }
    for (unsigned int i = 0; i < T.size(); i++)
    {
        this->chi2 += pow((this->Gdata.at(i) - this->CalcGfit(this->lambda, this->T.at(i))),2);
    }
    return;

}

double FreeEnergyFit::CalcGfit(vector <double> &lambda, double T)
{

    return lambda.at(0) + lambda.at(1)*(T-this->T0) + lambda.at(2)*T*log(T/this->T0);
}

double FreeEnergyFit::GetGfit(double T)
{
    return CalcGfit(this->lambda, T);
}

double FreeEnergyFit::GetTSfit(double T)
{
    return -T*( this->lambda[1] + this->lambda[2]*( 1 + log(T/this->T0) ) );
}

double FreeEnergyFit::GetHfit(double T)
{
    return GetGfit(T) + GetTSfit(T);
}

double FreeEnergyFit::ddchi2(double T, int i)
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

double FreeEnergyFit::calcf(vector <double> &lambda, int j)
{
    double sum = 0.0;
    for (unsigned int i = 0; i < this->T.size(); i++)
    {
        sum += (this->Gdata.at(i) - this->CalcGfit(lambda, this->T.at(i)) )*this->ddchi2(T.at(i), j);
    }
    return sum;
}

double FreeEnergyFit::GetLambda(int i)
{
    return this->lambda[i];
}

double FreeEnergyFit::GetLambdaGuess(int i)
{
    return this->lambda_init[i];
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

double FreeEnergyFit::GetT0()
{
    return this->T0;
}
