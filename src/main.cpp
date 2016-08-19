/*
 * (C) 2016 James W. Barnett
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "FreeEnergyFit.h"

using namespace std;

const double R = 8.314459848e-3; // kJ / K
const double CtoK = 273.15;

int main(int argc, char *argv[]) 
{
    if (argc != 2) 
    {
        cout << "Usage: " << endl;
        cout << "  " << argv[0] << " configfile" << endl;
        return -1;
    }
    const string configfile = argv[1];

    cout << "Reading " << configfile << "...";

    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini(configfile, pt);
    char *endptr;

    const double stepsize = strtod(pt.get<std::string>("fit.stepsize").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'fit.stepsize' must be a double." << endl;
    const double tol = strtod(pt.get<std::string>("fit.tol","1e-6").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'fit.tol' must be a double." << endl;
    const double T0_C = strtod(pt.get<std::string>("data.refT","298.15").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'data.refT' must be a double." << endl;
    const double T0 = T0_C + CtoK;
    const double alpha_init = strtod(pt.get<std::string>("fit.alpha").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'fit.alpha' must be a double." << endl;
    const double beta_init = strtod(pt.get<std::string>("fit.beta").c_str(), &endptr); 
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'fit.beta' must be a double." << endl;
    const double gamma_init = strtod(pt.get<std::string>("fit.gamma").c_str(), &endptr);
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'fit.gamma' must be a double." << endl;
    const int maxiter = strtol(pt.get<std::string>("fit.maxiter").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'fit.maxiter' must be an integer." << endl;
    const string outfile = pt.get<std::string>("output.file");
    const string units_in = pt.get<std::string>("units.in","kJ");
    const string units_out = pt.get<std::string>("units.out","kJ");
    const int temps_n = strtol(pt.get<std::string>("data.n").c_str(), &endptr, 10);
    if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'data.n' must be an integer." << endl;
    vector <double> T(temps_n);
    vector <double> T_C(temps_n);
    vector <double> G(temps_n);
    for (int i = 0; i < temps_n; i++)
    {
        T_C[i] = strtod(pt.get<std::string>("data.T"+to_string(i+1)).c_str(), &endptr);
        if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'data.T" + to_string(i+1) + "' must be a double." << endl;
        T[i] = T_C[i] + CtoK;
        G[i] = strtod(pt.get<std::string>("data.G"+to_string(i+1)).c_str(), &endptr);
        if (*endptr != ' ' && *endptr != 0) cout << "ERROR: 'data.G" + to_string (i+1) + "' must be a double." << endl;
    }
    cout << "done." << endl;

    double conv;
    if (units_in == units_out)
    {
        conv = 1.0;
    }
    else if (units_in == "kcal" && units_out == "kJ")
    {
        conv = 4.184;
    }
    else if (units_in == "kJ" && units_out == "kcal")
    {
        conv = 1.0/4.184;
    }
    else
    {
        cout << "ERROR: Units not recognized. Options are 'kcal' or 'kJ'. " << endl;
        return -1;
    }

    vector <double> x(3);
    x[0] = alpha_init;
    x[1] = beta_init;
    x[2] = gamma_init;
    FreeEnergyFit fit(x, T0, T, G, maxiter, stepsize, tol);
    fit.DoFit();

    ofstream oFS;
    oFS << scientific << setprecision(6);
    cout << "Writing output to " << outfile << "...";
    oFS.open(outfile.c_str());
    oFS << "# Free energy fit and its derivatives" << endl;
    oFS << "# Output file:          " << outfile << endl;
    oFS << "# T0 (K):               " << fit.GetT0() << endl;
    oFS << "# Max iterations:       " << fit.GetMaxiter() << endl;
    oFS << "# Step size:            " << fit.GetStepsize() << endl;
    oFS << "# Tolerance:            " << fit.GetTolerance() << endl;
    if (fit.IsConverged())
    {
        oFS << "# Fit converged" << endl;
    }
    else
    {
        oFS << "# WARNING: Fit did not converge." << endl;
    }
    oFS << "# Initial guesses: " << endl;
    oFS << "#   alpha   = " << fit.GetLambdaGuess(0) << endl;
    oFS << "#   beta    = " << fit.GetLambdaGuess(1) << endl;
    oFS << "#   gamma   = " << fit.GetLambdaGuess(2) << endl;
    oFS << "# Fitted parameters: " << endl;
    oFS << "#   alpha   = " << fit.GetLambda(0) << endl;
    oFS << "#   beta    = " << fit.GetLambda(1) << endl;
    oFS << "#   gamma   = " << fit.GetLambda(2)<< endl;
    oFS << "#   chi2    = " << fit.GetChi2() << endl;
    oFS << "#";
    oFS << setw(20) << "T (K) ";
    oFS << setw(20) << "Gdata (" + units_out + "/mol) ";
    oFS << setw(20) << "Gfit (" + units_out + "/mol) ";
    oFS << setw(20) << "-TS (" + units_out + "/mol) ";
    oFS << setw(20) << "H (" + units_out + "/mol) ";
    oFS << endl;
    for (int i = 0; i < temps_n; i++) 
    {
        oFS << setw(20) << fit.GetT(i);
        oFS << setw(20) << fit.GetGdata(i)*conv;
        oFS << setw(20) << fit.GetGfit(T[i])*conv;
        oFS << setw(20) << -fit.GetTSfit(T[i])*conv;
        oFS << setw(20) << fit.GetHfit(T[i])*conv;
        oFS << endl;
    }
    oFS.close();
    cout << "done." << endl;

    return 0;
}



