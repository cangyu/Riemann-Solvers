#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <exception>
#include <algorithm>

using namespace std;

const double a = 1.0;
const double ALPHA = 1.0;
const double BETA = 8.0;

const double CFL = 0.8;

double u0(double x)
{
    if (abs(x) <= 1.0)
        return ALPHA * exp(-BETA * pow(x, 2));
    else
        return 0.f;
}

const int NumOfPnt = 101;
const double xL = -1.0, xR = 1.0;
const double dx = (xR - xL)/(NumOfPnt - 1);
vector<double> x(NumOfPnt, 0.f);
vector<double> u_prev(NumOfPnt+2, 0.f);
vector<double> u_cur(NumOfPnt+2, 0.f);

const int NumOfStep = 12500;
const double dt = CFL * dx / a;

void output(ofstream &f, const vector<double> &x)
{
    const int n = x.size()-1;

    f << x[1];
    for (int i = 2; i < n; i++)
        f << '\t' << x[i];
    f << endl;
}

int main(int argc, char *argv[])
{
    cout << "===================================================" << endl;
    cout << "Numerical solution of the Linear Advection Problem:" <<endl; 
    cout << "\t u_t + a * u_x = 0"<<endl; 
    cout << "with smooth initial velocity profile."<<endl;

    //Grid coordinate & initial velocity profile
    x[0] = xL;
    for(int i = 1; i < NumOfPnt; i++)
        x[i] = x[i-1] + dx;

    //Scheme Coefs
    double b_l2 = 0, b_l1 = 0, b_0 = 0, b_r1 = 0;

    ofstream fout;

    cout << "========================CIR========================" << endl;
    fout.open("CIR.txt");
    if(!fout)
        throw "Failed to create output file!\n";

    fout << NumOfStep << "\t" << NumOfPnt << endl;
    for(int k = 0; k < NumOfPnt; k++)
        fout << x[k] << '\t';
    fout << endl;

    //IC
    for(int i = 1; i <= NumOfPnt; i++)
        u_prev[i] = u0(x[i-1]);
    
    //Periodical BC
    u_prev[0] = u_prev[NumOfPnt];
    u_prev[NumOfPnt+1] = u_prev[1];

    output(fout, u_prev);

    double a_plus = (a + abs(a))/2, a_minus = (a - abs(a))/2;
    double c_plus = dt * a_plus / dx, c_minus = dt * a_minus / dx;

    b_l1 = c_plus;
    b_0 = 1-abs(CFL);
    b_r1 = -c_minus;

    //Iterate over time
    for(int k = 1; k < NumOfStep; k++)
    {
        for(int i = 1; i <= NumOfPnt; i++)
            u_cur[i] = b_l1 * u_prev[i-1] + b_0 * u_prev[i] + b_r1 * u_prev[i+1];
        
        //Periodical BC
        u_cur[0] = u_cur[NumOfPnt];
        u_cur[NumOfPnt+1] = u_cur[1];

        output(fout, u_cur);
        
        u_prev.swap(u_cur);
    }
    fout.close();
    cout << "Done!" << endl;
    
    cout << "===================Lax-Friedrichs==================" << endl;
    fout.open("Lax-Friedrichs.txt");
    if(!fout)
        throw "Failed to create output file!\n";

    fout << NumOfStep << "\t" << NumOfPnt << endl;
    for(int k = 0; k < NumOfPnt; k++)
        fout << x[k] << '\t';
    fout << endl;

    //IC
    for(int i = 1; i <= NumOfPnt; i++)
        u_prev[i] = u0(x[i-1]);
    
    //Periodical BC
    u_prev[0] = u_prev[NumOfPnt];
    u_prev[NumOfPnt+1] = u_prev[1];

    output(fout, u_prev);

    b_l1 = (1+CFL)/2;
    b_0 = 0;
    b_r1 = (1-CFL)/2;

    //Iterate over time
    for(int k = 1; k < NumOfStep; k++)
    {
        for(int i = 1; i <= NumOfPnt; i++)
            u_cur[i] = b_l1 * u_prev[i-1] + b_r1 * u_prev[i+1];
        
        //Periodical BC
        u_cur[0] = u_cur[NumOfPnt];
        u_cur[NumOfPnt+1] = u_cur[1];

        output(fout, u_cur);
        
        u_prev.swap(u_cur);
    }
    fout.close();
    cout << "Done!" << endl;

    cout << "====================Lax-Wendroff===================" << endl;
    fout.open("Lax-Wendroff.txt");
    if(!fout)
        throw "Failed to create output file!\n";

    fout << NumOfStep << "\t" << NumOfPnt << endl;
    for(int k = 0; k < NumOfPnt; k++)
        fout << x[k] << '\t';
    fout << endl;

    //IC
    for(int i = 1; i <= NumOfPnt; i++)
        u_prev[i] = u0(x[i-1]);
    
    //Periodical BC
    u_prev[0] = u_prev[NumOfPnt];
    u_prev[NumOfPnt+1] = u_prev[1];

    output(fout, u_prev);

    b_l1 = CFL*(1+CFL)/2;
    b_0 = 1-pow(CFL, 2);
    b_r1 = -CFL*(1-CFL)/2;

    //Iterate over time
    for(int k = 1; k < NumOfStep; k++)
    {
        for(int i = 1; i <= NumOfPnt; i++)
            u_cur[i] = b_l1 * u_prev[i-1] + b_0 * u_prev[i] + b_r1 * u_prev[i+1];
        
        //Periodical BC
        u_cur[0] = u_cur[NumOfPnt];
        u_cur[NumOfPnt+1] = u_cur[1];

        output(fout, u_cur);
        
        u_prev.swap(u_cur);
    }
    fout.close();
    cout << "Done!" << endl;

    cout << "====================Warming-Beam===================" << endl;
    fout.open("Warming-Beam.txt");
    if(!fout)
        throw "Failed to create output file!\n";

    fout << NumOfStep << "\t" << NumOfPnt << endl;
    for(int k = 0; k < NumOfPnt; k++)
        fout << x[k] << '\t';
    fout << endl;

    //IC
    for(int i = 1; i <= NumOfPnt; i++)
        u_prev[i] = u0(x[i-1]);
    
    //Periodical BC
    u_prev[0] = u_prev[NumOfPnt];
    u_prev[NumOfPnt+1] = u_prev[1];

    output(fout, u_prev);

    b_l2 = CFL * (CFL - 1) / 2;
    b_l1 = CFL * (2 - CFL);
    b_0 = (CFL - 1) * (CFL - 2) /2;

    //Iterate over time
    for(int k = 1; k < NumOfStep; k++)
    {
        for(int i = 1; i <= NumOfPnt; i++)
            u_cur[i] = b_l2 * u_prev[i-2] + b_l1 * u_prev[i-1] + b_0 * u_prev[i];
        
        //Periodical BC
        u_cur[0] = u_cur[NumOfPnt];
        u_cur[NumOfPnt+1] = u_cur[1];

        output(fout, u_cur);
        
        u_prev.swap(u_cur);
    }
    fout.close();
    cout << "Done!" << endl;
    cout << "=========================End=======================" << endl;

    return 0;
}
