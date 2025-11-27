/// Use the Rk4 solver for coupled ODEs to solve for projectile 
/// motion with air resistance
///
/// Definition of our variables
/// x    = time <br>
/// y[0] = position along i axis  ; f_ri = dri/dt => velocity along i axis  <br>
/// y[1] = velocity along i axis  ; f_vi = dvi/dt => acceleration along i axis <br>
/// y[2] = position along j axis  ; f_rj = drj/dt => velocity along j axis <br>
/// y[3] = velocity along j axis  ; f_vj = dvj/dt => acceleration along j axis <br>


#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TText.h"

#include <functional>
#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;
namespace pl = std::placeholders;

constexpr double PI = 3.14159265;

struct Params {
  double g;   ///< acceleration [m/s^2]
  double m;   ///< mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double air_k;  ///< constant for air resistance. mass DOES matter with air resistance
} ;

// functions to describe simple projectile motion
// here use use ri,rj,rk to define directions to prevent confusion with
// standard ODE notation, where x=independent variable, \vec y=dependent variable(s)


/// \brief Change in position along \f$\hat i\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[1];
}

/// \brief Change in velocity along  \f$\hat i\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vi(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[1] / p->m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}

/// \brief Change in position along \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Air resistance model: F= \f$k v^2\f$
///
double f_rj(double x, const vector<double> &y, void *params=0){  
  (void) x;   // prevent unused variable warning
  return y[3];
}

/// Change in velocity along  \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vj(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[3] / p->m - p->g;
  // return -g;    // if no air constant acceleration along -j direction: F/m = -g
}

/// Total energy
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_e(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  return p->m*p->g*y[2] + 0.5*p->m*(y[1]*y[1] + y[3]*y[3]);
}

/// \brief Stopping condition
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Returns 0(1) to flag continuation(termination) of calculation 
double f_stop(double x, const vector<double> &y, void *params=0){
  (void) x;
  if (y[2]<0) return 1;  // stop calulation if the current step takes height to negative value
  return 0;  // continue calculation
}

/// \observer function for energy
TGraph *tgEnergy;
double f_eObs(double x, const vector<double> &y, void *params=0){
  tgEnergy->SetPoint(tgEnergy->GetN(), x, f_e(x, y, params));
  if (y[2]<0) return 1;  // stop calulation if the current step takes height to negative value
  return 0;  // continue calculation
}


//===========================================================================================================
//===========================================================================================================
//===========================================================================================================
//===========================================================================================================
//===========================================================================================================


/// \brief Use RK4 method to describe simple projectile motion.
int main(int argc, char **argv){

  // setup default parameters
  Params parsNoAir { 9.81, 10.0, 0.0 };
  void *p_parNoAir = (void*) &parsNoAir;

  double theta=PI/4;   // initial angle radians
  double v0=100;     // m/s

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************

  // *** test 2: Use RK4SolveN to calculate simple projectile motion
  vector<pfunc_t> v_fun{ f_ri, f_vi, f_rj, f_vj };   // 4 element vector of function pointers
  double xmax=20;     // tmax
  int nsteps=200;

  //Solve no air problem and find energy vs time
  tgEnergy = new TGraph();
  double xNoAir = 0;
  vector<double> yNoAir { 0, v0*cos(theta), 0, v0*sin(theta) };
  auto tgN_NoAir = RK4SolveN(v_fun, yNoAir, nsteps, xNoAir, xmax, p_parNoAir, f_eObs);

  int nRuns = 50;
  double massVec[nRuns];
  for(int i = 0; i <= nRuns; i++) 
    massVec[i] = ((10 - 0.5) / nRuns)*i + 0.5; //mass between 0.5kg and 10kg

  TGraph *tg_vterm = new TGraph();
  for(const auto m : massVec)
  {
    Params pars { 9.81, m, 0.1 };
    void *p_par = (void*) &pars;
    double x = 0;
    vector<double> y { 0, v0*cos(theta), 0, v0*sin(theta) };
    auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);
    tg_vterm->SetPoint(tg_vterm->GetN(), m, sqrt(y[1]*y[1] + y[3]*y[3]));
  }
  
  TCanvas *cNoAir = new TCanvas();
  tgN_NoAir[2].SetTitle("No air resistance rj vs. time;t (s);rj (m)");
  tgN_NoAir[2].Draw("al*");
  cNoAir->Draw();

  TCanvas *cNAEnergy = new TCanvas();
  tgEnergy->SetTitle("No air resistance energy vs. time;t (s); energy (J)");
  tgEnergy->Draw("al*");
  cNAEnergy->Draw();

  TCanvas *cvterm = new TCanvas();
  tg_vterm->SetTitle("Terminal velocity in air vs. Mass; mass (kg); vt (m/s)");
  tg_vterm->Draw("al*");
  cvterm->Draw();

  std::vector<std::string> comment {
      "The accuracy of the conservation of energy is improved with smaller step sized",
      "(i.e. more steps). Energy is conserved quite well at even n=200 steps, with only",
      "minor fluctuations in it, as can be seen in the graph.",
      "The solution for terminal velocity vs. mass when accounting for air resistance seem",
      "to be reasonably accurate. We notice both that the terminal velocity is always less than",
      "that of the no air resistance case and that it increases with larger masses",
      "- as we would expect by our physical intuition."
      };

  TCanvas* canvas = new TCanvas("c", "vterm Results", 800, 600);
  cNoAir->Print("vterm.pdf(");
  cNAEnergy->Print("vterm.pdf");
  cvterm->Print("vterm.pdf");
  float y = 0.9;
  for (const auto& line : comment) {
    TText* t = new TText(0.1, y, line.c_str());
    t->SetTextSize(0.03);
    t->Draw();
    y -= 0.06; // spacing between lines
  }
  canvas->Print("vterm.pdf)");

  cout << "Terminal velocity with no air (m=10kg) = " << sqrt(yNoAir[1]*yNoAir[1]+yNoAir[3]*yNoAir[3]) << " m/s" << endl;
  cout << "Terminal velocity with air (m=10kg) = " << tg_vterm->GetY()[39]  << " m/s" << endl;

  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

