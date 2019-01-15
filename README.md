# Riemann Solvers
Code snippets follow from ___Riemann Solvers and Numerical Methods for Fluid Dynamics___ by Eleuterio F. Toro, where essentials of CFD are discussed in detail.

## Linear Advection(ch2 & ch5)
Target equation:
$$u_t + a u_x = 0$$  
Both _smooth_ and _discontinous_ initial velocity profile are examined. The exact solution is quite trival, just trace back along the line $\frac{dx}{dt} = a$ so that $\frac{du}{dt} = 0$. 
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
Only the __discontinous__ initial velocity profile is examined, and analytically, the exact solution is either a _shock_ wave($U_L > U_R$) or a _rarefaction_ wave($U_L \le U_R$).
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

## Euler Equation(theory in ch3)
1-D Euler equation:
$$U_t + F(U)_x = 0$$

where
$$U = [\rho, \rho u, E]^T$$

and
$$F = [\rho u, \rho u^2 + p, u(E + p)]^T$$

($E = \rho(\frac{1}{2}u^2 + e)$, the total energy per volumn)  

For ideal gases, the specific internal energy $e$ is expressed as
$$e = \frac{p}{(\gamma-1)\rho}$$

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

