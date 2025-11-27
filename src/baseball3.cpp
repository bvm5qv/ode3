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

// constants
constexpr double PI = 3.14159265;
constexpr double g = 9.81;
constexpr double B = 4.1E-4;

// Parameters
struct Params {
  double g;   // acceleration [m/s^2]
  double d;   // m diameter of ball
  double B;   // B param for air resistance
  double w;   // magnitude of angular velocity
  double phi;   // angle of angular velocity vector (rad)
};

// data vectors to hold results
std::vector<double> vX{};
std::vector<double> vY{};
std::vector<double> vZ{};
std::vector<double> vdX{};
std::vector<double> vdY{};
std::vector<double> vdZ{};

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

//drag force equation (in m/s)
double F_drag(double v) { return 0.0039+0.0058/(1+std::exp((v-35)/5)); } 

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
  return -F_drag(v)*v*S[vardX] + p->B*p->w*(S[vardZ]*sin(p->phi) - S[vardY]*cos(p->phi));
}
//acceleration in y direction
double f_dy(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  Params *p = (Params*)params;
  double v = sqrt(S[vardX]*S[vardX] + S[vardY]*S[vardY] + S[vardZ]*S[vardZ]);
  return -F_drag(v)*v*S[vardY] + p->B*p->w*S[vardX]*cos(p->phi);
}
//acceleration in z direction
double f_dz(double t, const vector<double> &S, void *params=0){ 
  (void) t;   // prevent unused variable warning
  Params *p = (Params*)params;
  double v = sqrt(S[vardX]*S[vardX] + S[vardY]*S[vardY] + S[vardZ]*S[vardZ]);
  return -g - F_drag(v)*v*S[vardZ] - p->B*p->w*S[vardX]*sin(p->phi);
}

//===========================================================================================
//===========================================================================================
//===========================================================================================


/// Returns 0(1) to flag continuation(termination) of calculation and records data
double xend;
double f_obs(double t, const vector<double> &S, void *params=0){
  (void) t;
  vX.push_back(S[varX]);
  vY.push_back(S[varY]);
  vZ.push_back(S[varZ]);
  vdX.push_back(S[vardX]);
  vdY.push_back(S[vardY]);
  vdZ.push_back(S[vardZ]);
  //std::cout << "(x,y,z) = " << "(" << S[varX] << ", " << S[varY] << ", " << S[varZ] << ")\n";
  if (S[varX]>xend) return 1;  // stop calulation after passing xmax
  return 0;  // continue calculation
}

void SetupSlider(vector<double>& y0, Params *p)
{
  double v0 = 38;
  double theta0 = PI/180;
  y0 = { 0, 0, 0, v0*cos(theta0), 0, v0*sin(theta0) };
  p->g = g;
  p->d = 0.01;   // m
  p->B = B;
  p->w = 30*2*PI;     // s^{-1}
  p->phi = 0;    // rad
}

void SetupCurve(vector<double>& y0, Params *p)
{
  double v0 = 38;
  double theta0 = PI/180;
  y0 = { 0, 0, 0, v0*cos(theta0), 0, v0*sin(theta0) };
  p->g = g;
  p->d = 0.01;    // m
  p->B = B;
  p->w = 30*2*PI;      // s^{-1}
  p->phi = PI/4;  // rad
}

void SetupScrewball(vector<double>& y0, Params *p)
{
  double v0 = 38;
  double theta0 = PI/180;
  y0 = { 0, 0, 0, v0*cos(theta0), 0, v0*sin(theta0) };
  p->g = g;
  p->d = 0.01;      // m
  p->B = B;
  p->w = 30*2*PI;        // s^{-1}
  p->phi = 3*PI/4;  // rad
}

void SetupFastball(vector<double>& y0, Params *p)
{
  double v0 = 42.5;
  double theta0 = PI/180;
  y0 = { 0, 0, 0, v0*cos(theta0), 0, v0*sin(theta0) };
  p->g = g;
  p->d = 0.01;      // m
  p->B = B;
  p->w = 30*2*PI;        // s^{-1}
  p->phi = 5*PI/4;  // rad
}

struct Trajectory {
    std::vector<double> x, y, z, dx, dy, dz;
    TGraph *gY;
    TGraph *gZ;
};

Trajectory ComputeTrajectory(void (*SetupFn)(vector<double>&, Params*)) {
    vector<double> y0;
    Params *p = new Params();
    SetupFn(y0, p);

    // reset global trackers
    vX.clear(); vY.clear(); vZ.clear();
    vdX.clear(); vdY.clear(); vdZ.clear();

    double t = 0;
    double tmax = 60;
    int nsteps = 2000;
    vector<pfunc_t> v_fun{ f_x, f_y, f_z, f_dx, f_dy, f_dz };

    xend = 60 / 3.3;
    auto tg = RK4SolveN(v_fun, y0, nsteps, t, tmax, p, f_obs);

    int N = vX.size();
    Trajectory out;
    out.x.resize(N); out.y.resize(N); out.z.resize(N);
    out.dx.resize(N); out.dy.resize(N); out.dz.resize(N);

    for (int i = 0; i < N; i++) {
        out.x[i] = 3.28084 * vX[i];
        out.y[i] = 3.28084 * vY[i];
        out.z[i] = 3.28084 * vZ[i];
        out.dx[i] = 3.28084 * vdX[i];
        out.dy[i] = 3.28084 * vdY[i];
        out.dz[i] = 3.28084 * vdZ[i];
    }

    out.gY = new TGraph(N, out.x.data(), out.y.data());
    out.gZ = new TGraph(N, out.x.data(), out.z.data());

    return out;
}

int main(int argc, char **argv) {

    TApplication theApp("theApp", &argc, argv);

    // List of pitch types
    std::vector<std::pair<std::string, void(*)(std::vector<double>&, Params*)>> pitchSetups = {
        {"Slider",     SetupSlider},
        {"Curveball",  SetupCurve},
        {"Screwball",  SetupScrewball},
        {"Fastball",   SetupFastball}
    };

    TCanvas *c = new TCanvas("c", "All pitches", 1200, 800);
    c->Divide(2, 2);

    int pad = 1;
    for (auto &p : pitchSetups) {
        auto traj = ComputeTrajectory(p.second);

        c->cd(pad++);
        traj.gY->SetTitle((p.first + " Trajectory; x (ft); y/z (ft)").c_str());
        traj.gY->SetLineColor(kBlue);
        traj.gZ->SetLineColor(kRed);

        traj.gY->SetMinimum(-4);
        traj.gY->SetMaximum(2);
        traj.gY->Draw("AL");
        traj.gZ->Draw("L SAME");

        auto leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(traj.gY, "y", "l");
        leg->AddEntry(traj.gZ, "z", "l");
        leg->Draw();
    }

    c->SaveAs("pitches.pdf");

    std::cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
    return 0;
}


