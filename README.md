# Riemann Solvers
Code snippets follow from ___Riemann Solvers and Numerical Methods for Fluid Dynamics___ by Eleuterio F. Toro, where essentials of CFD are discussed in detail.

## Linear Advection
Different schemes are employed for comparison
> * CIR  
> * Lax-Friedrichs  
> * Lax-Wendroff  
> * Warming-Beam  
> * Godunov  

Usage:  
> * Compile: `g++ smooth.cc -o a.out` or `g++ discontinuous.cc -o a.out` 
> * Execute: `./a.out`  
> * Plot: `python animate.py`

## Invisid Burgers Equation
Different schemes are employed for comparison
> * CIR  
> * Lax-Friedrichs  
> * Lax-Wendroff  
> * Warming-Beam  
> * Godunov  

Usage:  
> * Compile: `g++ main.cc -o a.out`
> * Execute: `./a.out`  
> * Plot: `python animate_single.py` or `python animate_all.py`

## Euler Equation
1-D Euler equation with ideal gas.
### Exact solution
Calculate numerically the exact solution.  
Especially, vacuum condition is considered.
> * Compile: `g++ main.cc -o Exact.out`  
> * Execute: `./Exact.out < inp.dat`  
> * Plot: `python animate.py`
### Godunov Scheme
Calculate the approximate solution where the exact solution of Riemann problem is applied locally.  
> * Compile: `g++ main.cc -std=c++11 -o Godunov.out`  
> * Execute: `./Godunov.out < inp.dat`  
> * Plot: `python animate.py`

