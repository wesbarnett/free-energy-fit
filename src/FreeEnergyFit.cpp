
#include "FreeEnergyFit.h"

FreeEnergyFit::FreeEnergyFit(vector <double> &lambda, double T0, vector <double> &T, vector <double> &G, int max_iter, double stepsize, double tol)
{

    this->T0 = T0;
    this->xdata = T;
    this->ydata = G;
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
    return;
}

double FreeEnergyFit::CalcFit(vector <double> &lambda, double T)
{

    return lambda.at(0) + lambda.at(1)*(T-this->T0) + lambda.at(2)*T*log(T/this->T0);
}

double FreeEnergyFit::GetTSfit(double T)
{
    return -T*( this->lambda[1] + this->lambda[2]*( 1 + log(T/this->T0) ) );
}

double FreeEnergyFit::GetHfit(double T)
{
    return GetFit(T) + GetTSfit(T);
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

double FreeEnergyFit::GetT0()
{
    return this->T0;
}
