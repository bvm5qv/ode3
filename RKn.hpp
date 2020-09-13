// generic solvers for linear first order ODE
// for convenience these function return a graph of the solution for y(x)
#pragma once 
#include <TGraph.h>
#include <vector>
using std::vector;

// generic fcn pointer for a linear differential equation, f = y'(x,y)
// double (*f)(double x, double y)
// the equation may depend on one independent variable (x) and a
// dependent variable (y)
// y0: initial condition
// nsteps, x0, xmax are used to se t the range and h, the step size
TGraph RK4Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax);

/**********************************************
RK4 solver for set of N 1st order ODEs (1 independent parameter)
The list of eqns is stored in a vector of fcn pointers
initial conditions are given in a vector of doubles

B/c the 1st order ODEs may be coupled, each equation will have 
access to the entire set of dependent parameters, thus we use 
the interface: double func_t(double X, const vector<double> &Y),
where X represents the independent paramater and Y is a vector of
the dependent parameters 

The solver will return a set of graphs for each dependent variable
***********************************************/
typedef double func_t(double, const vector<double>&);
typedef func_t* pfunc_t;   // function pointer matching our ODE interface

// fnlist: vector of function pointers to the ODEs describing the system
// y0: vector of initial conditions (need one per ODE)
// nsteps: maximum number of steps in simulation
// x0,xmax: starting and maximum value of dependent vaiable for simulation
// fstop: optional function of dependent varibles to define some stopping
//        condition based on the results of the simulation
vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
			 int nsteps, double x0, double xmax, pfunc_t fstop=0);

// this version is intended to run until the stopping condition is met
// nmax: limit the number of steps taken in case stopping isn't satisfied 
vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
			 double h, double x0, pfunc_t fstop, int nmax=1000);


// this is a version of the algorithm with an adaptive step size
// errdef: is requested error per step, the error calculation
//         corresponds to a relative (absolute) error estimate
//         for y values large (small) compared to errdef
// maxrep: max recalulations of h at each step to prevent infinite loops
vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
			  int nsteps, double x0, double xmax, pfunc_t fstop=0,
			  double errdef=1e-9, int maxrep=5);

// this version is intended to run until the stopping condition is met
// fstop function is required!
// maxsteps: limit number of steps taken if stopping condition isn't satisfied 
vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
			  double h,  double x0, pfunc_t fstop,
			  double errdef=1e-9, int maxrep=5, int maxsteps=1000);

