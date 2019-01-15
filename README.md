# Riemann Solvers
Code snippets follow from ___Riemann Solvers and Numerical Methods for Fluid Dynamics___ by Eleuterio F. Toro, where essentials of CFD are discussed in detail.

## Linear Advection(ch2 & ch5)
Both _smooth_ and _discontinous_ initial velocity profile are examined. The exact solution is quite trival, just trace back along the characteristic line.  
Different schemes are employed for comparison:
> * CIR  
> * Lax-Friedrichs  
> * Lax-Wendroff  
> * Warming-Beam  
> * Godunov  

Usage:  
> * Compile: `g++ smooth.cc -o a.out` or `g++ discontinuous.cc -o a.out` 
> * Execute: `./a.out`  
> * Plot: `python3 animate.py`

## Invisid Burgers Equation(ch2 & ch5)
Target equation:
$$u_t + u u_x = 0$$  
Only the __discontinous__ initial velocity profile is examined, and analytically, the exact solution is either a _shock_ wave or a _rarefaction_ wave.
Different schemes are employed for comparison:
> * CIR  
> * Lax-Friedrichs  
> * Lax-Wendroff  
> * Warming-Beam  
> * Godunov  

Usage:  
> * Compile: `g++ main.cc -o a.out`
> * Execute: `./a.out`  
> * Plot: `python3 animate_single.py` or `python3 animate_all.py`

## Euler Equation
1-D Euler equation(ch3) with ideal gases.

### Exact solution(ch4)
The general solution of the exact solution follows the 3-wave pattern, where the _contact_ must lies in between, _shock_ and _rarefaction_ waves stay at left or right.  
The exact solution can be calculate numerically, where a iterative procedure is necessary for solving the _pressure_. The exact solution requies much computational effort and this is why approximate riemann solvers are studied extensively back in 1980s.
Especially, vacuum condition is considered for completeness.

Usage:  
> * Compile: `g++ main.cc -o Exact.out`  
> * Execute: `./Exact.out < inp.dat`  
> * Plot: `python3 animate.py`

### Godunov Scheme(ch6)
Calculate the approximate solution where the exact solution of Riemann problem is applied locally.  
> * Compile: `g++ main.cc -std=c++11 -o Godunov.out`  
> * Execute: `./Godunov.out < inp.dat`  
> * Plot: `python3 animate.py`

