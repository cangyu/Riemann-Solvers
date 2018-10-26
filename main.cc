#include <iostream>
#include <sstream>
#include <iomanip>
#include <exception>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

const double TOL = 1e-6;

const double gamma = 1.4;
const double G1 = 0.5*(gamma - 1) / gamma;
const double G2 = 0.5*(gamma + 1) / gamma;
const double G3 = 2 * gamma / (gamma - 1);
const double G4 = 2 / (gamma - 1);
const double G5 = 2 / (gamma + 1);
const double G6 = (gamma - 1) / (gamma + 1);
const double G7 = (gamma - 1) / 2;
const double G8 = gamma - 1;
const double G9 = -G2;
const double G10 = -0.5*(gamma + 1) / pow(gamma, 2);
const double G11 = -0.5*(3 * gamma + 1) / gamma;

typedef struct
{
	double rho;
	double u;
	double p;

	//Auxiliary variables introduced for efficiency and convinence
	double a;
}PrimitiveVariable;

inline double sound_speed(double p, double rho)
{
	return sqrt(gamma*p / rho);
}

inline double Ak(const PrimitiveVariable &Wk)
{
	return G5 / Wk.rho;
}

inline double Bk(const PrimitiveVariable &Wk)
{
	return G6 * Wk.p;
}

inline double gk(double p, const PrimitiveVariable &Wk)
{
	return sqrt(Ak(Wk) / (p + Bk(Wk)));
}

double fk(double p, const PrimitiveVariable &Wk)
{
	if (p > Wk.p)
	{
		//shock
		return (p - Wk.p)*gk(p, Wk);
	}
	else
	{
		//rarefraction
		return G4 * Wk.a*(pow(p / Wk.p, G1) - 1);
	}
}

double f(double p, const PrimitiveVariable &Wl, const PrimitiveVariable &Wr)
{
	return fk(p, Wl) + fk(p, Wr) + (Wr.u - Wl.u);
}

double dfk(double p, const PrimitiveVariable &Wk)
{
	if (p > Wk.p)
	{
		//shock
		double A = Ak(Wk);
		double B = Bk(Wk);
		return sqrt(A / (p + B))*(1 - 0.5*(p - Wk.p) / (B + p));
	}
	else
	{
		//rarefraction
		return pow(p / Wk.p, G9) / (Wk.rho*Wk.a);
	}
}

double df(double p, const PrimitiveVariable &Wl, const PrimitiveVariable &Wr)
{
	return dfk(p, Wl) + dfk(p, Wr);
}

double ddfk(double p, const PrimitiveVariable &Wk)
{
	if (p > Wk.p)
	{
		//shock
		double A = Ak(Wk);
		double B = Bk(Wk);
		return -0.25 * sqrt(A / (p + B)) * (4 * B + 3 * p + Wk.p) / pow(B + p, 2);
	}
	else
	{
		//rarefraction
		return G10 * Wk.a / pow(Wk.p, 2) * pow(p / Wk.p, G11);
	}
}

double ddf(double p, const PrimitiveVariable &Wl, const PrimitiveVariable &Wr)
{
	return ddfk(p, Wl) + ddfk(p, Wr);
}

int main(int argc, char *argv[])
{
	PrimitiveVariable Wl, Wr;

	int n = 0;
	cin >> n;
	for (int i = 1; i <= n; i++)
	{
		//Input parameters
		try
		{
			cin >> Wl.rho >> Wl.u >> Wl.p;
			cin >> Wr.rho >> Wr.u >> Wr.p;
		}
		catch (...)
		{
			cerr << "Invalid input!" << endl;
			return -1;
		}

		if (!(Wl.rho > 0) || !(Wr.rho > 0))
		{
			cerr << "Invalid data!" << endl;
			return -2;
		}

		Wl.a = sound_speed(Wl.p, Wl.rho);
		Wr.a = sound_speed(Wr.p, Wr.rho);

		stringstream ss;
		ss << i;
		string header = "Case " + ss.str();

		cout.fill('=');
		cout << setw(32) << std::right << header << setw(32) << std::right << "" << endl;
		cout.fill(' ');

		//Check pressure positivity
		double du = Wr.u - Wl.u;
		double du_crit = G4 * (Wl.a + Wr.a);
		if (!(du_crit > du))
			throw "Vacuum!";

		//Select initial pressure
		double p_min = min(Wl.p, Wr.p);
		double p_max = max(Wl.p, Wr.p);
		double f_min = f(p_min, Wl, Wr);
		double f_max = f(p_max, Wl, Wr);

		cout << setw(16) << std::left << "P_L" << setw(16) << std::left << "P_R" << setw(16) << std::left << "f_min" << setw(16) << std::left << "f_max" << endl;
		cout << setw(16) << std::left << Wl.p << setw(16) << std::left << Wr.p << setw(16) << std::left << f_min << setw(16) << std::left << f_max << endl;
		cout << endl;

		double p_m = (Wl.p + Wr.p) / 2;
		p_m = max(TOL, p_m);
		double p_pv = p_m - du * (Wl.rho + Wr.rho)*(Wl.a + Wr.a) / 8;
		p_pv = max(TOL, p_pv);
		double gL = gk(p_pv, Wl);
		double gR = gk(p_pv, Wr);
		double p_ts = (gL * Wl.p + gR * Wr.p - du) / (gL + gR);
		p_ts = max(TOL, p_ts);
		double p_tr = pow((Wl.a + Wr.a - G7 * du) / (Wl.a / pow(Wl.p, G1) + Wr.a / pow(Wr.p, G1)), G3);
		p_tr = max(TOL, p_tr);

		cout << setw(16) << std::left << "P_TR" << setw(16) << std::left << "P_PV" << setw(16) << std::left << "P_TS" << setw(16) << std::left << "P_M" << endl;
		cout << setw(16) << std::left << p_tr << setw(16) << std::left << p_pv << setw(16) << std::left << p_ts << setw(16) << std::left << p_m << endl;
		cout << endl;

		cout << "P0: ";
		double p0 = p_m;
		if (f_min < 0 && f_max < 0)
		{
			cout << "Two-Shock Approximation" << endl;
			p0 = p_ts;
		}
		else if (f_min > 0 && f_max > 0)
		{
			cout << "Two-Rarefraction Approximation" << endl;
			p0 = p_tr;
		}
		else
		{
			cout << "Intermediate Approximation" << endl;
			p0 = p_pv;
		}

		//Solve pressure
		cout << "\nSolving p* ..." << endl;
		int iter_cnt = 0;
		double CHA = 1.0;
		while (CHA > TOL)
		{
			++iter_cnt;

			double fder = df(p0, Wl, Wr);
			if (fder == 0)
				throw "Zero derivative!";

			double fval = f(p0, Wl, Wr);
			double fder2 = ddf(p0, Wl, Wr);
			double p = p0 - fval * fder / (pow(fder, 2) - 0.5*fval*fder2);
			if (p < 0)
			{
				p0 = TOL;
				break;
			}

			CHA = abs(2 * (p - p0) / (p + p0));
			p0 = p;
			cout << "Iter: " << iter_cnt << "  " << "P: " << p0 << endl;
		}

		cout << "Done!\n" << endl;

		double us = (Wl.u + Wr.u) / 2 + (fk(p0, Wr) - fk(p0, Wl))/2;
		double rhosL = 0.0;
		double rhosR = 0.0;

		cout << setw(16) << std::left << "p*" << setw(16) << std::left << "u*" << setw(16) << std::left << "rho_sL" << setw(16) << std::left << "rho_sR" << endl;
		cout << setw(16) << std::left << p0 << setw(16) << std::left << us << setw(16) << std::left << rhosL << setw(16) << std::left << rhosR << endl;
		cout << endl;
	}

	cout.fill('=');
	cout << setw(32) << std::right << "End" << setw(32) << std::right << "" << endl;
	cout.fill(' ');

	return 0;
}