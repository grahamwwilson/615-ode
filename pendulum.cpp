//
// pendulum.cpp 
//
//                   Graham W. Wilson   
//                       updated May 2018
//                       updated October 2019 for HW3 F2019
//
// Solve initial value problem ODE of finite angle simple pendulum
// 
//    d(theta)/dt = omega
//    d(omega)/dt = -g/l sin(theta)
// 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <string>
using namespace std;

enum algorithms { EULER, RK2, RK4 };
     
int main (int argc, char *argv[])
{
//
// Command line arguments are: 
// time-step, out-filename, dat-filename, 
// stepping algorithm.
//
// Stepping algorithm choices are:
// 0 = Euler
// 1 = Second order Runge Kutta (RK2)
// 2 = Fourth order Runge Kutta (RK4) 
// (they are defined by the algorithms enumeration)
//
   cout << "argc = " << argc << endl;
   for (int i=0; i<argc;++i){
       cout << "argv [ " << i << "] = " << argv[i] << endl;
   }
   double dt;
   if(argc > 1)dt = atof(argv[1]);
   cout << "dt = " << dt << endl;

   ofstream outfile;
   outfile.open(argv[2]);

   ofstream datfile;
   datfile.open(argv[3]);

   int algorithm = atoi(argv[4]);

   cout << "Algorithm choice set to " << algorithm << endl;


// Use reasonably precise value of pi computed trigonometrically
   double pie = 4.0*atan(1.0);
   cout << "pie = " << setprecision(16) << pie << endl;
  
// Problem constants
   const double g=9.81;
   const double l=1.5;
   const double mass = 1.0;   // Not really needed but define it anyway

// Initial conditions
   double theta = pie/3.0;
   double omega = 0.0;

   double ktheta[4],komega[4];     // RK4 corrections to theta, omega

   double energy = mass*l*(0.5*l*omega*omega + g*(1.0-cos(theta))); 
   double e0 = energy;   // Initial energy for reference

  //
  // Use FORTRAN compatibility output using examples from 
  // page 29 of Barton and Nackmann, Scientific and Engineering C++
  //
   outfile << setiosflags(ios::showpoint | ios::uppercase);
   outfile << setiosflags(ios::fixed);
   outfile.precision(16);

   datfile << setiosflags(ios::showpoint | ios::uppercase);
   datfile << setiosflags(ios::fixed);
   datfile.precision(16);

   if(algorithm == EULER){
      cout << "Using Euler stepping algorithm" << endl;
      outfile << "Using Euler Method with time step of " 
              << dt << endl;
   }
   else if(algorithm == RK2){
      cout << "Using RK2 stepping algorithm" << endl;
      outfile << "Using RK2 with time step of " 
              << dt << endl;
   }
   else if(algorithm == RK4){
      cout << "Using RK4 stepping algorithm" << endl;
      outfile << "Using 4th order Runge-Kutta Method with time step of " 
              << dt << endl;
   }
   else{
      cout << "No defined stepping algorithm...." << endl;
   }

   // Initialize time stepping
   int istep = 0;
   double t = 0.0;

   // Initial Values
   outfile << "Step"      << setw(8)  << istep 
           << " t: "      << setw(20) << t 
           << " theta: "  << setw(20) << theta 
           << " omega: "  << setw(20) << omega 
           << " energy: " << setw(20) << energy <<  endl;

   datfile << setw(20) << t 
           << setw(20) << theta 
           << setw(20) << omega 
           << setw(20) << 1.0 - (energy/e0) <<  endl;

   while(t<10.0){ 

      ktheta[0] = dt*omega;
      komega[0] = dt*( (-g/l)*sin(theta) ) ;

      ktheta[1] = dt*(omega + 0.5*komega[0]);
      komega[1] = dt*( (-g/l)*sin(theta + 0.5*ktheta[0]) ) ;

      ktheta[2] = dt*(omega + 0.5*komega[1]);
      komega[2] = dt*( (-g/l)*sin(theta + 0.5*ktheta[1]) ) ;

      ktheta[3] = dt*(omega + komega[2]);
      komega[3] = dt*( (-g/l)*sin(theta + ktheta[2]) ) ;

    // Increment t and number of steps
      istep++;
      t = double(istep)*dt;

    // Implement the appropriate correction formula for both variables
      if (algorithm == EULER){
      // Euler
         theta = theta + ktheta[0];
         omega = omega + komega[0];
      }
      else if (algorithm == RK2){
      // 2nd order Runge-Kutta
         theta = theta + ktheta[1];
         omega = omega + komega[1]; 
      }
      else if (algorithm == RK4){
      // 4th order Runge-Kutta
         theta = theta + ( ktheta[0] + 2.0*ktheta[1] + 2.0*ktheta[2] + ktheta[3] )/ 6.0;
         omega = omega + ( komega[0] + 2.0*komega[1] + 2.0*komega[2] + komega[3] )/ 6.0; 
      }
      else{
         cout << "Still no defined stepping algorithm...." << endl;
      }

      // Re-evaluate energy
      energy = mass*l*(0.5*l*omega*omega + g*(1.0-cos(theta)));    
      double efrac = 100.0*(energy-e0)/e0;
      if(istep%1000==0){
         cout << "Step " << istep << " Percentage Deviation in Energy " 
              << setw(20) << efrac << " (%)" << endl;
      }

      outfile << "Step"      << setw(8)  << istep 
              << " t: "      << setw(20) << t 
              << " theta: "  << setw(20) << theta 
              << " omega: "  << setw(20) << omega 
              << " fenergy: " << setw(20) << (energy/e0)-1.0 <<  endl;

      datfile << setw(20) << t 
              << setw(20) << theta 
              << setw(20) << omega 
              << setw(20) << 1.0 - (energy/e0) <<  endl;
   }
   outfile.close();
   datfile.close();

   cout << " " << endl;

   return 0;
}
