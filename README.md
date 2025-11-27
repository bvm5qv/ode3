# ode3

#Atticus Yohn

physx630 odelib
---

To build the ODE library and example programs, simply type `make` in this top level ode3 directory.

Description of example programs:<br>

**RKnTest**: Solves a single 1st order ODE using the single equation RK4 solver and the ODE array solver

**RKnStep**: A basic example of the ODE array solver is applied to projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>. At each step in the elapsed time and x,y positions are printed.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m


**RKnDemo**: Solves for projectile motion with a simple model of air resistance, force of air resistance = -kv^2<br>
This program includes graphical output.  Detailed output is saved in TGraph objects in RKnDemo.root.  The file **RKnPlotDemo.py** shows how to access date in the TGraphs and can be used to generate additional plots.<br>
Optional parameters [default values]<br>
* -v initial_velocity [100] m/s
* -t angle_thera [45] degrees
* -m mass_of_projectile [10] kg
* -k coefficient_of_air_resistance [0.1] kg/m

**baseball1**:  Starter template for first baseball problem

**baseball2**:  Starter template for second baseball problem

**baseball_drag.ipynb**: this notebook describes the drag force equations used in the text.

gsl starter code
---

The starter code here (projGSL.cpp) demonstrates very basic usage of the gsl for solving a problem of coupled differential equations. An 8th order R-K solver with fixed step size is used. You are encouraged to try other solvers as you explore the problem. See here for the gsl docs: https://www.gnu.org/software/gsl/doc/html/ode-initval.html

This example solves the 2D projectile motion problem with a simple model for air resistance. After each step, data are stored in ROOT TGraphs, which are then displayed at the conclusion of the calculation.

The gsl provides a number of ODE solvers and a variety of interfaces.  Some of the solvers (not R-K methods) use the Jacobian matrix, which gives the devivative of the function wrt the dependent parameters.  See the gsl examples for details.

Python starter code
---

Two examples are given for using ODE solvers from the scipy.integrate sub-package in Python. In these examples graphs are made using matplotlib.

    Solution (projScPY2.py[ipynb]) using a more modern interface scipy.integrate.solve_ivp. See also: https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html and https://www.programcreek.com/python/example/119375/scipy.integrate.solve_ivp

    Solution (projScPY.py[ipynb]) using an older interface scipy.integrate.odeint¶ (see comments here: https://docs.scipy.org/doc/scipy/reference/integrate.html).  I do not recommend using this interface any longer.

The notebook versions contain additional comments on using the integrators.

# =================================================================================================
# =================================================================================================
# ===========================================  Answers ============================================
# =================================================================================================
# =================================================================================================


RESULT OF QUESTION 1 PART B: (baseball1.cpp output)
diameter d = 0.007500 m
mass m = 0.145000 kg
x(0) = 0 m
y(0) = 0 m
z(0) = 1.400000 m
vx(0) = 45.155122 m/s
vy(0) = 0.000000 m/s
vz(0) = 0.788186 m/s
theta_0 = 1.000000 degrees
xend(tmax) = 18.517503 m
yend(tmax) = 0.000000 m
zend(tmax) = 0.897356 m
********************************
(xend,z0,theta0) = (18.500000,1.400000,1.000000)
v_pitch = 45.162000 m/s
********************************


NOTES ON QUESTION 2: (baseball2.cpp and baseball3.cpp)
In this problem, I edited the starter code (baseball2.cpp) and reached the solution without changing the overall framework of the code or ending lines. The code as there written returns one single trajectory graph relating to one single type of pitch (depending on the value of ip which must be changed by hand). Given the sloppy way in which I wrote the code, I was not sure how to compute all cases and display them in a single pdf by running the code a single time. Consequently, I refactored the code in baseball3.cpp, which is capable of running all four cases for the purpose of creating the desired graphs (though it destroys the end data and is more disimilar from the format of the starter code). In any event, both versions of the code are technically "fully functioning", I just created and used baseball3.cpp to more easily generate the pitches.pdf image. 


