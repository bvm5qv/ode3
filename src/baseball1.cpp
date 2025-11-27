///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

constexpr double PI = 3.14159265;
constexpr double g = 9.81;
constexpr double b = 1.6E-4;
constexpr double c = 0.25;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};

enum : unsigned char //index of variables in our dependant variable vector (S)
{
  varX,    // pos X
  varY,    // pos Y
  varZ,    // pos Z
  vardX,   // dXdt
  vardY,   // dYdt
  vardZ,   // dZdt

  TOTAL_NUM_VARS
};

//===========================================================================================
//====================== Differential Equations =============================================
//===========================================================================================

//rate of change in x direction
double f_x(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  return S[vardX];
}
//rate of change in y direction
double f_y(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  return S[vardY];
}
//rate of change in z direction
double f_z(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  return S[vardZ];
}

//acceleration in x direction
double f_dx(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  Params *p = (Params*)params;
  double v = sqrt(S[vardX]*S[vardX] + S[vardY]*S[vardY] + S[vardZ]*S[vardZ]);
  double Fnet = p->b*p->d*v + p->c*p->d*p->d*v*v; // magnitude of force
  return (Fnet / p->m) * (-S[vardX] / v);
}
//acceleration in y direction
double f_dy(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  Params *p = (Params*)params;
  double v = sqrt(S[vardX]*S[vardX] + S[vardY]*S[vardY] + S[vardZ]*S[vardZ]);
  double Fnet = p->b*p->d*v + p->c*p->d*p->d*v*v; // magnitude of force
  return (Fnet / p->m) * (-S[vardY] / v);
}
//acceleration in z direction
double f_dz(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  Params *p = (Params*)params;
  double v = sqrt(S[vardX]*S[vardX] + S[vardY]*S[vardY] + S[vardZ]*S[vardZ]);
  double Fnet = p->b*p->d*v + p->c*p->d*p->d*v*v; // magnitude of force
  return (Fnet / p->m) * (-S[vardZ] / v) - p->g;
}

//===========================================================================================
//===========================================================================================
//===========================================================================================

/// \brief Stopping condition
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Returns 0(1) to flag continuation(termination) of calculation 
double f_stop(double t, const vector<double> &S, void *params=0){
  (void) t;
  //std::cout << "(x,y,z) = " << "(" << S[varX] << ", " << S[varY] << ", " << S[varZ] << ")\n";
  if (S[varZ]<0.9) return 1;  // stop calulation if the current step takes height below strike zone
  return 0;  // continue calculation
}

int main(int argc, char **argv){

  // examples of parameters
  Params pars { g, 0.145, 0.0075, b, c };
  void *p_par = (void*) &pars;

  double xend = 18.5;  // meters to plate
  double z0 = 1.4;     // height of release [m]
  double theta0 = 1;   // angle of velocity at release (degrees)
                                      
  bool showPlot=false;    // keep this flag false by default
  
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  double vPitch = 0;   // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here

  double vPitchStep = 0.001;
  double xFinal = 0;
  int nsteps = 1000;
  vector<pfunc_t> v_fun{ f_x, f_y, f_z, f_dx, f_dy, f_dz };
  vector<double> y;
  while( xend - xFinal > 0.01 )
  { 
    vPitch += vPitchStep;
    double t = 0;
    double tmax = 2*xend / vPitch;
    vector<double> y0 { 0, 0, z0, vPitch*cos(theta0 * PI/180), 0, vPitch*sin(theta0 * PI/180) };
    auto tg = RK4SolveN(v_fun, y0, nsteps, t, tmax, p_par, f_stop);
    xFinal = y0[varX];
    y = y0;
    //std::cout << "vPitch = " << vPitch << "\nxFinal = " << xFinal << "\nzFinal = " << y0[varZ] << "\n\n";
  }

  printf("diameter d = %lf m\n", pars.d);
  printf("mass m = %lf kg\n", pars.m);
  printf("x(0) = 0 m\n");
  printf("y(0) = 0 m\n");
  printf("z(0) = %lf m\n", z0);
  printf("vx(0) = %lf m/s\n", vPitch*cos(theta0 * PI/180)); 
  printf("vy(0) = %lf m/s\n", 0.0); 
  printf("vz(0) = %lf m/s\n", vPitch*sin(theta0 * PI/180)); 
  printf("theta_0 = %lf degrees\n", theta0); 
  printf("xend(tmax) = %lf m\n", y[varX]);
  printf("yend(tmax) = %lf m\n", y[varY]);
  printf("zend(tmax) = %lf m\n", y[varZ]);

  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

