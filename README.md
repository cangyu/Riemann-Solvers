# Riemann Solvers
Code snippets follow from ___Riemann Solvers and Numerical Methods for Fluid Dynamics___ by Eleuterio F. Toro, where essentials of CFD are discussed in detail.

## Linear Advection(ch2 & ch5)
Both _smooth_ and _discontinous_ initial velocity profile are examined.  
The exact solution is quite trival, just trace back along the characteristic line.  
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
Only the __discontinous__ initial velocity profile is examined.  
Analytically, the exact solution is either a _shock_ wave or a _rarefaction_ wave.
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

### Godunov's Method(ch6)
The essential ingredient of Godunov's method is to solve _Riemann Problem_ locally, and the keypoint in numerical parctice is to identify all of the 10 possible wave patterns so that the inter-cell flux can be calculated properly.

Usage:
> * Compile: `g++ main.cc -std=c++11 -o Godunov.out`  
> * Execute: `./Godunov.out < inp.dat`  
> * Plot: `python3 animate.py`

### Lax-Friedrichs(ch5 & ch6)
Follow the guide to rewrite the differential scheme into conservative form.  
It's interesting to notice that Lax-Friedrichs scheme is identical to the Riemann Solution averaged at the __half__ of each time step.
The key part in parctice is the Lax-Friedrichs inter-cell flux, see (5.77).

Usage:
> * Compile: `g++ main.cc -std=c++11 -o Lax.out`  
> * Execute: `./Lax.out < inp.dat`  
> * Plot: `python3 animate.py`

### Richtmyer(ch5 & ch6)
Again, rewrite the differential scheme into its conservative form.  
Lax-Wendroff is of 2nd-order accuracy in space and time.  
The Richtmyer version(2-step Lax-Wendroff), introduces an intermediate step where a tempory conservative variable is computed, and the inter-cell flux is calculated later on based on this tempory conservative variable.

Usage:
> * Compile: `g++ main.cc -std=c++11 -o Richtmyer.out`  
> * Execute: `./Richtmyer.out < inp.dat`  
> * Plot: `python3 animate.py`

### RCM(ch7)
TODO

### FVS(ch8)
Here, the inter-cell flux is not calculatd directly from the exact solution of Riemann Problem. Instead, the flux at each point is splitted into 2 parts: __upstream traveling__ and __downstream traveling__, then, for each inter-cell, the flux is seen as the sum of the upstream traveling part from the __left__ point and the downstream traveling part from the __right__ point.  
3 typical splitting techniques are introduced.  

#### Steger-Warming

Usage:
> * Compile: `g++ main.cc -std=c++11 -o SW.out`  
> * Execute: `./SW.out < inp.dat`  
> * Plot: `python3 animate.py`

#### van Leer

Usage:
> * Compile: `g++ main.cc -std=c++11 -o vL.out`  
> * Execute: `./vL.out < inp.dat`  
> * Plot: `python3 animate.py`

#### AUSM

Usage:
> * Compile: `g++ main.cc -std=c++11 -o AUSM.out`  
> * Execute: `./AUSM.out < inp.dat`  
> * Plot: `python3 animate.py`

### Approximate-state(ch9)
Godunov's method in conjunction with exact Riemann solution spends lots of computational effort on solving pressure in star region iteratively. One kind of approximate Riemann solution seeks to approximate the __state__ of the wave patterns so that unkonwn values at star region can be solved immediately.  
5 typical approximate-state Riemann Solvers are implemented.  

#### PVRS

Usage:
> * Compile: `g++ main.cc -std=c++11 -o PVRS.out`  
> * Execute: `./PVRS.out < inp.dat`  
> * Plot: `python3 animate.py`

#### TRRS

Usage:
> * Compile: `g++ main.cc -std=c++11 -o TRRS.out`  
> * Execute: `./TRRS.out < inp.dat`  
> * Plot: `python3 animate.py`

#### TSRS

Usage:
> * Compile: `g++ main.cc -std=c++11 -o TSRS.out`  
> * Execute: `./TSRS.out < inp.dat`  
> * Plot: `python3 animate.py`

#### AIRS

Usage:
> * Compile: `g++ main.cc -std=c++11 -o AIRS.out`  
> * Execute: `./AIRS.out < inp.dat`  
> * Plot: `python3 animate.py`

#### ANRS

Usage:
> * Compile: `g++ main.cc -std=c++11 -o ANRS.out`  
> * Execute: `./ANRS.out < inp.dat`  
> * Plot: `python3 animate.py`

### HLL(ch10)
By esitmate the 2 fastest wave spreading speed, the averaged flow variable in between are __constant__!

#### Direct estimation under Roe-average
The wave spreading speed at left and right front are esitmated according to corresponding __eigenvalues__. The key part in practice is the __Roe-averaged eigenvalues__.  

Usage:
> * Compile: `g++ direct.cc -std=c++11 -o Euler.out`  
> * Execute: `./Euler.out < inp.dat`  
> * Plot: `python3 animate.py`

#### Pressure-based estimation
Approximate the pressure at star region, and then estimate the speed at wave front from the exact solution of Riemann Problem.  

Usage:
> * Compile: `g++ pressure_based.cc -std=c++11 -o Euler.out`  
> * Execute: `./Euler.out < inp.dat`  
> * Plot: `python3 animate.py`

#### Rusanov
Only the maximun wave spreading speed is approximated according to eigenvalues, this is a much simpler choice for the fastest signal velocities.  

Usage:
> * Compile: `g++ rusanov.cc -std=c++11 -o Euler.out`  
> * Execute: `./Euler.out < inp.dat`  
> * Plot: `python3 animate.py`

### HLLC(ch10)
An internal wave is added compared with HLL so that contact is better resolved.

#### Direct estimation under Roe-average
Same as that in HLL.  

Usage:
> * Compile: `g++ direct.cc -std=c++11 -o Euler.out`  
> * Execute: `./Euler.out < inp.dat`  
> * Plot: `python3 animate.py`

#### Pressure-based estimation
Same as that in HLL.  

Usage:
> * Compile: `g++ pressure_based.cc -std=c++11 -o Euler.out`  
> * Execute: `./Euler.out < inp.dat`  
> * Plot: `python3 animate.py`

### Roe(ch11)
Instead of estimating the signal prompting speed, another approach seeks to approximate the Jacobian matrix with known left and right states such that 3 essential properties(Hyperbolicity, Consistency and Conservation across discontinuities) are satisfied. An ingenious way to construct such a approximated Jacobian matrix is given by Roe, where the famous parameter vector consists the __square root of density__ is introduced. In practice, the inter-cell flux is computed from the wave strength, eigenvalues and eigenvectors explicitly.

Usage:
> * Compile: `g++ main.cc -std=c++11 -o Euler.out`  
> * Execute: `./Euler.out < inp.dat`  
> * Plot: `python3 animate.py`

### Osher(ch12)
TODO
