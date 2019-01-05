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

inline double kinetic_energy(double u)
{
    return 0.5 * pow(u, 2);
}

class PrimitiveVar
{
private:
    double rho, u, p;
    double a, e;

public:
    PrimitiveVar(double density, double velocity, double pressure)
    {
        set(density, velocity, pressure);
    }

    PrimitiveVar(istream &in)
    {
        double density, velocity, pressure;
        in >> density >> velocity >> pressure;
        set(density, velocity, pressure);
    }

    PrimitiveVar() 
    {
        set(1.0, 0.0, 101325.0);
    }

    ~PrimitiveVar() {}

    void set(double density, double velocity, double pressure)
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

//Exact solution of then 1-D Euler equation
//This is the Riemann problem, where initial discontinuity exists
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

inline double E(const PrimitiveVar &W)
{
    return W.rho * (internal_energy(W.p, W.rho) + kinetic_energy(W.u));
}

inline double E(double density, double velocity, double pressure)
{
    return density * (internal_energy(pressure, density) + kinetic_energy(velocity));
}

class ConservativeVar
{
private:
    double u[3];

public:
    ConservativeVar() 
    {
        set(1.0, 0.0, 101325.0);
    }

    ConservativeVar(const PrimitiveVar &x)
    {
        set(x.rho, x.u, x.p);
    }

    ~ConservativeVar() {}

    void set(double density, double velocity, double pressure)
    {
        u[0] = density;
        u[1] = density * velocity;
        u[2] = E(density, velocity, pressure);
    }
};

class FluxVar
{
private:
    double f[3];

public:
    FluxVar() 
    {
        f[0] = f[1] = f[2] = 0.0;
    }

    ~FluxVar() {}

    void set(const PrimitiveVar &W)
    {
        f[0] = W.rho * W.u;
        f[1] = W.rho * pow(W.u, 2) + W.p;
        f[2] = W.u * (E(W) + W.p);
    }

	void set(double density, double velocity, double pressure)
	{
		f[0] = density * velocity;
		f[1] = density * pow(velocity, 2) + pressure;
		f[2] = velocity * (E(density, velocity, pressure) + W.p);
	}

	double operator[](int n)
	{
		if(n < 0 || n >=3)
			throw("Invalid index!");
		return f[n];
	}
};

void W_fanL(PrimitiveVar *W, double S, PrimitiveVar &ans)
{
    double density = W->rho * pow(G5 + G6 / W->a *(W->u - S), G4);
	double velocity = G5 * (W->a + G7 * W->u + S);
	double pressure = W->p * pow(G5 + G6 / W->a * (W->u - S), G3);
    ans.set_param(density, velocity, pressure);
}

void W_fanR(PrimitiveVar *W, double S, PrimitiveVar &ans)
{
	double density = W->rho * pow(G5 - G6 / W->a *(W->u - S), G4);
	double velocity = G5 * (-W->a + G7 * W->u + S);
	double pressure = W->p * pow(G5 - G6 / W->a * (W->u - S), G3);
    ans.set_param(density, velocity, pressure);
}

class InterCell
{
private:
	PrimitiveVar *Wl, *Wr;
	double p_s, u_s, rho_sL, rho_sR;
	FluxVar local_flux;

public:
	InterCell()
	{
		Wl = NULL;
		Wr = NULL;
		p_s = 101325.0;
		u_s = 0.0;
		rho_sL = 1.0;
		rho_sR = 1.0;
	}

	~InterCell() {}

	void set_param(PrimitiveVar *left, PrimitiveVar *right)
	{
		Wl = left;
		Wr = right;

		//Get the exact solution
		p_s = p_star(*left, *right);
		u_s = u_star(p_s, *left, *right);
		rho_sL = rho_star(p_s, *left);
		rho_sR = rho_star(p_s, *right);

		PrimitiveVar x;
		RP(0.0, x);
		local_flux.set_param(x);
	}

    //Solution of the Riemann Problem
    //10 possible wave patterns
    void RP(double S, PrimitiveVar &ans)
    {
        if(S < u_s)
        {
            //At the left of the contact discontinuity
            if(p_s > Wl->p)
            {
                //Left shock
                double S_L = Wl->u - Wl->a * sqrt(G2 * p_s / Wl->p + G1);
                if(S < S_L)
                    ans = *Wl;
                else
                    ans.set_param(rho_sL, u_s, p_s);
            }
            else
            {
                //Left fan
                double S_HL = Wl->u - Wl->a;
                if (S < S_HL)
                    ans = *Wl;
                else
                {
                    double S_TL = u_s - sound_speed(p_s, rho_sL);
                    if(S > S_TL)
                        ans.set_param(rho_sL, u_s, p_s);
                    else
                        W_fanL(Wl, S, ans);
                }
            }
        }
        else
        {
            //At the right of the contact discontinuity
            if(p_s > Wr->p)
            {
                //Right shock
                double S_R = Wr->u + Wr->a * sqrt(G2 * p_s / Wr->p + G1);
                if(S < S_R)
                    ans.set_param(rho_sR, u_s, p_s);
                else
                    ans = *Wr;
            }
            else
            {
                //Right fan
                double S_HR = Wr->u + Wr->a;
                if(S > S_HR)
                    ans = *Wr;
                else
                {
                    double S_TR = u_s + sound_speed(p_s, rho_sR);
                    if(S < S_TR)
                        ans.set_param(rho_sR, u_s, p_s);
                    else
                        W_fanR(Wr, S, ans);
                }
            }
        }
    }
};

const int NumOfPnt = 101;
const double xL = 0, xR = 1.0;
const double xM = (xL + xR) / 2;
const double dx = (xR - xL) / (NumOfPnt - 1);

void output(ofstream &f, int time_step, const vector<PrimitiveVar> &W)
{
	f << time_step << endl;
	for (int i = 0; i < NumOfPnt; i++)
		f << W[i].rho << '\t' << W[i].u << '\t' << W[i].p << endl;
}

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
		fout << x[0];
		for (int i = 1; i < NumOfPnt; i++)
			fout << '\t' << x[i];
		fout << endl;


		//Initialize
		int PREV = 0, CUR = 1;
		vector<vector<PrimitiveVar>> w(2, vector<PrimitiveVar>(NumOfPnt+2));
		w[PREV][0] = Wl;
		for(int k = 1; k <= NumOfPnt; ++k)
		{
			if(x[k] < xM)
				w[PREV][k] = Wl;
			else
				w[PREV][k] = Wr;
		}
		w[PREV][NumOfPnt+1] = Wr;
		output(fout, 0, W[PREV]);

		vector<InterCell> f(NumOfPnt+1);

		//Iterate
		for (int k = 1; k < NumOfStep; ++k)
		{
			for(int j = 0; j < NunOfPnt+1; ++j)
				f[j].set_param(&w[k], &w[j+1]);

			for(int j = 1; j <=NumOfPnt; ++j)
			{
				double factor = - dt/dx;
				double fc[3];
				for(int i = 0; i < 3; i++)
					fc[i] = f[j].local_flux[i] - f[j-1].local_flux[i];
				
				w[CUR][j].density = w[PREV][j].density + factor * fc[0];
				w[CUR][j].velocity = w[PREV][j].velocity + factor * fc[1];
				w[CUR][j].pressure = w[PREV][j].pressure + factor * fc[2];
			}
			output(fout, k, w[CUR]);
		}
	}
	
	return 0;
}
