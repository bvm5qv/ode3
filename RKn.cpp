#include "RKn.hpp"
#include <assert.h>
#include <algorithm>
#include <cmath>

using std::max;
using std::min;

TGraph RK4Solve(double (*f)(double x, double y), double y0,
		int nsteps, double x0, double xmax){

  double h=(xmax-x0)/nsteps;     // step size
  double x=x0;                   // independent variable
  
  double y=y0;          // dependent variable to plot vs x
  TGraph tg;
  tg.SetPoint(0,x0,y0);
	      
  double k1,k2,k3,k4;
  for (int i=0; i<nsteps-1; i++){
    k1 = h * f(x,y);
    k2 = h * f(x+h/2,y+k1/2);
    k3 = h * f(x+h/2,y+k2/2);
    k4 = h * f(x+h,y+k3);
    y = y + k1/6 + k2/3 + k3/3 + k4/6;
    x+=h;
    tg.SetPoint(i+1,x,y);
  }
  return tg;
}

// called by RK4SolveN, this function performs one step using the RK4 algorithm
// and returns the new vector of y values
vector<double> RK4StepN(vector<pfunc_t> &f, vector<double> y,
			double x, double h){
  int nFcn=f.size();
  vector<double> k1(nFcn), k2(nFcn), k3(nFcn), k4(nFcn);
  vector<double> ytmp(nFcn);

  for (int i=0; i<nFcn; i++){
    k1[i] = h * f[i](x, y);
    ytmp[i] = y[i] + k1[i]/2;
  }
  for (int i=0; i<nFcn; i++){
    k2[i] = h * f[i](x+h/2, ytmp);
    ytmp[i] = y[i] + k2[i]/2;
  }
  for (int i=0; i<nFcn; i++){
    k3[i] = h * f[i](x+h/2, ytmp);
    ytmp[i] = y[i] + k3[i];
  }
  for (int i=0; i<nFcn; i++){
    k4[i] = h * f[i](x+h, ytmp);
  }
  // calculate next step
  for (int i=0; i<nFcn; i++) {
    y[i] = y[i] + k1[i]/6 + k2[i]/3 + k3[i]/3 + k4[i]/6;
  }
  return y;
}

// make the graphs that will return information on the solution
vector<TGraph> SetupGraphs(int nFcn){
  vector<TGraph> tg(nFcn);
  for (int i=0; i<nFcn; i++){
    char buf[100];
    sprintf (buf,";independent variable (x);dependent variable y[%d]",i);
    tg[i].SetTitle(buf);
    sprintf(buf,"xy%d",i);
    tg[i].SetName(buf);
  }
  return tg;
}


vector<TGraph> RK4SolveN(vector<pfunc_t> &f, vector<double> &y0,
			 int nsteps, double x0, double xmax, pfunc_t fstop){

  assert(f.size() == y0.size());   // need one initial condition per function
  int nFcn=f.size();
  double h=(xmax-x0)/nsteps;     // step size
  double x=x0;                   // independent variable
  vector<double> y = y0;         // start with initial conditions on dependent vars
  vector<double> ytmp(nFcn);
  
  vector<TGraph> tg = SetupGraphs(nFcn);
  for (int i=0; i<nFcn; i++) tg[i].SetPoint(0,x0,y0[i]);
  
  for (int n=0; n<nsteps-1; n++){
    ytmp=RK4StepN(f, y, x, h);
    if (fstop && fstop(x+h,ytmp)) break;
    // advance to next step and store results in graphs
    y=ytmp;
    x+=h;
    for (int i=0; i<nFcn; i++) tg[i].SetPoint(n+1,x,y[i]);
  }
  
  return tg;
}


vector<TGraph> RK4SolveN(vector<pfunc_t> &fnlist, vector<double> &y0,
			 double h, double x0, pfunc_t fstop,
			 int nmax){
  double xmax=x0+h*nmax;  // serves as a backup stopping condition 
  return RK4SolveN(fnlist, y0, nmax, x0, xmax, fstop);
}


// this code loosely follows the examples in Fitzpatrick and
// Numerical Recipes to adapt the step sizes
vector<TGraph> RK4SolveNA(vector<pfunc_t> &f, vector<double> &y0,
			  int nsteps, double x0, double xmax, pfunc_t fstop,
			  double errdef, int maxrep){
  // **
  // The following could be passed as parameters, to allow more control
  // when calling the function but here they are hardcoded to (hopefully)
  // resonable values to simplify the function interface
  const double hmin=errdef*4;
  const double hmax=(xmax-x0)/4;
  // **

  assert(f.size() == y0.size());   // need one initial condition per function
  int nFcn=f.size();
  double h=(xmax-x0)/nsteps;     // step size
  double h_last=0, h_new=0;      // for tracking adapted step size
  double x=x0;                   // independent variable
  vector<double> y = y0;         // start w/ init. conditions on dependent vars
  vector<double> y1,y2;          // temporary storage for adapting steps
  
  vector<TGraph> tg = SetupGraphs(nFcn);
  for (int i=0; i<nFcn; i++) tg[i].SetPoint(0,x0,y0[i]);

  while (x<xmax){ // loop over steps in the x range
    int nreps=0;
    
    while (nreps<maxrep){ // take the steps, adapting step size as needed
      //full step
      y1=RK4StepN(f, y, x, h);
      // two half steps
      y2=RK4StepN(f, y, x, h/2);
      y2=RK4StepN(f, y2, x+h/2, h/2);
      h_last=h;     // last step size used to evaluate the functions
      
      // Calculate truncation error, comparing the full step to two half steps
      double errAbs=0, errRel=0;
      for (int i=0; i<nFcn; i++){
	errAbs = max(errAbs,fabs(y1[i] - y2[i]));
	errRel = max(errRel,fabs(errAbs / y2[i]));
      }
      double err=max(errAbs,errRel);  // error estimate
      err=max(err,1e-15);            // protect against unlikely case of 0 error
      h_new = h * pow (fabs (errdef / err), 0.2);
      h_new=max(h_new,hmin);    // restrict adapted step size to above limits
      h_new=min(h_new,hmax);
      if (err<errdef) {
	h=h_new;
	break;                  // this step satisfies the requested errdef
      }
      nreps++;
      h=h_new;                  // if not try again
    }
    
    
    if (fstop && fstop(x+h_last,y2)) break;  // fix me
    // advance to next step and store results in graphs
    y=y2;
    for (int i=0; i<nFcn; i++) tg[i].SetPoint(tg[i].GetN(),x+h_last,y[i]);
    x+=h_new;
  }

  return tg;
}

vector<TGraph> RK4SolveNA(vector<pfunc_t> &fnlist, vector<double> &y0,
			  double h, double x0, pfunc_t fstop,
			  double errdef, int maxrep, int maxsteps){
  double xmax=x0+h*maxsteps;  // serves as a backup stopping condition
  return RK4SolveNA(fnlist, y0, maxsteps, x0, xmax, fstop, errdef, maxrep);
}
