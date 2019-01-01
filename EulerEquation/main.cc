#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>

const double G0 = 1.4;
const double G1 = 0.5*(G0 - 1) / G0;
const double G2 = 0.5*(G0 + 1) / G0;
const double G3 = 2 * G0 / (G0 - 1);
const double G4 = 2 / (G0 - 1);
const double G5 = 2 / (G0 + 1);
const double G6 = (G0 - 1) / (G0 + 1);
const double G7 = (G0 - 1) / 2;
const double G8 = G0 - 1;
const double G9 = -G2;
const double G10 = -0.5*(G0 + 1) / pow(G0, 2);
const double G11 = -0.5*(3 * G0 + 1) / G0;
const double G12 = 1 / G0;

using namespace std;

inline double sound_speed(double p, double rho)
{
    return sqrt(G0 * p / rho);
}

inline double internal_energy(double p, double rho)
{
    return p / (G8 * rho);
}

class PrimitiveVar
{
public:
    double rho, u, p;
    double a, e;

public:
    PrimitiveVar(double density, double velocity, double pressure)
    {
        set_param(density, velocity, pressure);
    }

    PrimitiveVar(istream &in)
    {
        double density, velocity, pressure;
        in >> density >> velocity >> pressure;
        set_param(density, velocity, pressure);
    }

    PrimitiveVar() 
    {
        set_param(1.0, 0.0, 101325.0);
    }

    ~PrimitiveVar() {}

    void set_param(double density, double velocity, double pressure)
    {
        rho = density;
        u = velocity;
        p = pressure;
        a = sound_speed(p, rho);
        e = internal_energy(p, rho);
    }
};

inline double Ak(const PrimitiveVar &W)
{
	return G5 / W.rho;
}

inline double Bk(const PrimitiveVar &W)
{
	return G6 * W.p;
}

inline double gk(double p, const PrimitiveVar &W)
{
	return sqrt(Ak(W) / (p + Bk(W)));
}

inline double fk(double p, const PrimitiveVar &W)
{
    if(p > W.p)
        return (p - W.p) * gk(p, W);
    else
        return G4 * W.a * (pow(p/W.p, G1) - 1);
}

inline double dfk(double p, const PrimitiveVar &W)
{
    if(p > W.p)
    {
        double B = Bk(W);
        return sqrt(Ak(W)/(B + p)) * (1 - 0.5 * (p - W.p) / (B + p));
    }
    else
        return 1.0 / (W.rho * W.a) / pow(p/W.p, G2);
}

inline double ddfk(double p, const PrimitiveVar &W)
{
    if(p > W.p)
    {
        double B = Bk(W);
        return -0.25*sqrt(Ak(W)/(B+p))*(4*B+3*p+W.p)/pow(B+p,2);
    }
    else
        return G10*W.a/pow(W.p, 2)*pow(p/W.p, G11);
}

inline double f(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
    return fk(p, Wl) + fk(p, Wr) + (Wr.u - Wl.u);
}

inline double df(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
    return dfk(p, Wl) + dfk(p, Wr);
}

inline double ddf(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
    return ddfk(p, Wl) + ddfk(p, Wr);
}

double p_star(const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
    const double TOL = 1e-6;
    const double du = Wr.u - Wl.u;

    //Select initial pressure
    const double f_min = f(min(Wl.p, Wr.p), Wl, Wr);
    const double f_max = f(max(Wl.p, Wr.p), Wl, Wr);

    double p_m = (Wl.p + Wr.p) / 2;
    p_m = max(TOL, p_m);

    double p_pv = p_m - du * (Wl.rho + Wr.rho)*(Wl.a + Wr.a) / 8;
    p_pv = max(TOL, p_pv);

    const double gL = gk(p_pv, Wl);
    const double gR = gk(p_pv, Wr);
    double p_ts = (gL * Wl.p + gR * Wr.p - du) / (gL + gR);
    p_ts = max(TOL, p_ts);

    double p_tr = pow((Wl.a + Wr.a - G7 * du) / (Wl.a / pow(Wl.p, G1) + Wr.a / pow(Wr.p, G1)), G3);
    p_tr = max(TOL, p_tr);

    double p0 = p_m;
    if (f_min < 0 && f_max < 0)
        p0 = p_ts;
    else if (f_min > 0 && f_max > 0)
        p0 = p_tr;
    else
        p0 = p_pv;

    //Solve
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
    }

    return p0;
}

inline double u_star(double p, const PrimitiveVar &Wl, const PrimitiveVar &Wr)
{
    return (Wl.u + Wr.u) / 2 + (fk(p, Wr) - fk(p, Wl)) / 2;
}

inline double rho_star(double p, const PrimitiveVar &W)
{
	double t = p / W.p;

	if (t > 1.0)
		return W.rho * ((t + G6) / (G6*t + 1));
	else
		return W.rho * pow(t, G12);
}

class InterCellPnt
{
private:
    PrimitiveVar Wl, Wr;
    double p_s, u_s, rho_sL, rho_sR;

public:
    InterCellPnt(const PrimitiveVar &left, const PrimitiveVar &right):Wl(left), Wr(right)
    {
        //Get the exact solution
        p_s = p_star(left, right);
        u_s = u_star(p_s, left, right);
        rho_sL = rho_star(p_s, left);
        rho_sR = rho_star(p_s, right);
    }

    ~InterCellPnt() {}
};

const int NumOfPnt = 101;
const double xL = 0, xR = 1.0;
const double xM = (xL + xR) / 2;
const double dx = (xR - xL) / (NumOfPnt - 1);

int main(int argc, char **argv)
{
	int n;
    double dt;
    int NumOfStep;

    //Coordinates
    vector<double> x(NumOfPnt, xL);
    for(int k = 1; k < NumOfPnt; ++k)
        x[k] = x[k-1] + dx;

    //Loop all cases
	cin >> n;
	for(int k = 0; k < n; ++k)
	{    
        //Input
		PrimitiveVar Wl(cin), Wr(cin);
        cin >> dt >> NumOfStep;
        
        //Output coordinates and intial settings
		ofstream fout("godunov.txt");
		if (!fout)
			throw "Failed to open file!";

		fout << NumOfStep << '\t' << NumOfPnt << endl;
		for (int i = 0; i < NumOfPnt; i++)
			fout << x[i] << endl;
		fout << 0 << endl;
		for (int i = 0; i < NumOfPnt; i++)
		{
			if (x[i] < xM)
				fout << Wl.rho << '\t' << Wl.u << '\t' << Wl.p << endl;
			else
				fout << Wr.rho << '\t' << Wr.u << '\t' << Wr.p << endl;
		}

        //Solve
        //TODO
	}
	
	return 0;
}