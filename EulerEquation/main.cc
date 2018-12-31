#include <fstream>
#include <cmath>

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

double sound_speed(double p, double rho)
{
    return sqrt(G0 * p / rho);
}

double internal_energy(double p, double rho)
{
    return p / (G8 * rho);
}

class PrimitiveVar
{
private:
    double rho, u, p;
    double a, e;

public:
    PrimitiveVar(double density, double velocity, double pressure)
    {
        rho = density;
        u = velocity;
        p = pressure;
        a = sound_speed(p, rho);
        e = internal_energy(p, rho);
    }

    ~PrimitiveVar() {}
};

void RiemannSolver(const PrimitiveVar &left, const PrimitiveVar &Wr, double S, double &p_s, double &u_s, double &rho_sL, double &rho_sR)
{

}


class InterCellPnt
{
private:
    PrimitiveVar Wl, Wr;
    double p_s, u_s, rho_sL, rho_sR;

public:
    InterCellPnt(const PrimitiveVar &left, const PrimitiveVar &right):Wl(left), Wr(right)
    {
        RiemannSolver(left, right, 0, p_s, u_s, rho_sL, rho_sR);
    }

    ~InterCellPnt() {}
};

int main(int argc, char **argv)
{
    ifstream fin('inp.txt');

    fin.close();

    return 0;
}