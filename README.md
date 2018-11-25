# RiemannSolver
Riemann Solvers and Numerical Methods for Fluid Dynamics 

## Exact Solver
Solve the 1D Euler equation with ideal gas.  
Especially, vacuum condition is considered.
> * Compile: `g++ main.cc -o Exact.out`  
> * Execute: `./Exact.out < inp.txt`  
> * Plot: `python animate.py`

## Linear Advection
Different schemes are employed for comparison
> * CIR  
> * Lax-Friedrichs  
> * Lax-Wendroff  
> * Warming-Beam  
> * Godunov  

Usage:  
> * Compile: `g++ smooth.cc -o a.out`  or `g++ discontinuous.cc -o a.out` 
> * Execute: `./a.out`  
> * Plot: `python animate.py`
