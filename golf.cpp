#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
using namespace std;
using namespace boost::numeric::odeint;

ofstream out, dout;
const double G = 9.81;  // g in m s^-2 
//const double RHO = 1.225; // kg m^-3 at 101325 Pa (1 atm)
//const double RHO = 1.1773; // kg m^-3 at 101325 Pa (1 atm) at 80F  
const double RHO = 1.18; // kg m^-3
const double R = 0.84*2.54e-2; // m
const double MASS = 1.62/35.274; // kg

// Mere mortal
const double V0 = 150.0*0.44704; // in m/s
const double THETA0 = 12.0*M_PI/180.0;  // in rad
const double SPINRATE = 2750.0; // this is in RPM

// Tiger Woods
//const double V0 = 180.0*0.44704; // in m/s
//const double THETA0 = 11.0*M_PI/180.0;  // 11 degrees
//const double SPINRATE = 2200.0; // this is in RPM

typedef boost::array< double , 4 > state_type;

// The rhs of x' = f(x,t)
// For [x] = {x, y, vx, vy}
void rhs( const state_type &x , state_type &dxdt , const double t )
{

//these ones are easy
    dxdt[0] =  x[2];                // dx/dt  = vx
    dxdt[1] =  x[3];                // dy/dt  = vy

    double vsq = x[2]*x[2] + x[3]*x[3];
    double S = (SPINRATE*2.0*M_PI/60.0)*R/sqrt(vsq);
    double area = M_PI*R*R;
    
//    double cD = 0.1995028 + 0.1889649*S + 1.4650386*S*S;
    double cD = 0.199 + 0.189*S + 1.47*S*S;
 
    double kdrag = 0.5*cD*RHO*area;

    double Fdrag = kdrag*vsq;
    double Fdragx = -Fdrag*x[2]/sqrt(vsq);
    double Fdragy = -Fdrag*x[3]/sqrt(vsq);

//    double cL = -3.25*S*S + 1.99*S;   // Fig 4b from Lyu, Kensrud, Smith and Tosaya.
//    double cL = 0.0693859 + 0.9879173*S;
    double cL = 0.0694 + 0.988*S;

    double FMagnus = 0.5*cL*RHO*area*vsq;
    double FMagnusx = - FMagnus*x[3]/sqrt(vsq);  //  
    double FMagnusy =   FMagnus*x[2]/sqrt(vsq);  // If x is horizontal, lift is upwards
    
    dxdt[2] = (Fdragx + FMagnusx)/MASS;
    dxdt[3] = -G + (Fdragy + FMagnusy)/MASS;
 
}

void write_out( const state_type &x , const double t )
{
    double ypm = 1.0/0.9144;
    out  << setw(5) << t << " " << x[0] << " " << x[1] << " " << x[2] 
              << " " << x[3] << endl;

    double v = sqrt(x[2]*x[2] + x[3]*x[3]);
    v = v/0.44704; // in mph
    dout  << setw(5) << t << " " << x[0]*ypm << " " << x[1]*ypm << " " 
          << x[2]*ypm << " " << x[3]*ypm << " " << v << endl; 
}

typedef runge_kutta_fehlberg78< state_type > stepper_type;

int main()
{
    state_type x = { 0.0, 0.0, V0*cos(THETA0), V0*sin(THETA0) }; // initial conditions

    out.open("golf.out");
    out << "Time(s)    x(m)   y(m)   vx(m/s)   vy(m/s) " << endl; 
    out.precision(16);
    dout.open("golf.dat");
    dout << "Time (s)   x(yards)  y(yards)  vx(yards/s)  vy(yards/s)  v(mph) " << endl;
    dout.precision(16);
    integrate_n_steps( make_controlled(1e-14,1e-14,stepper_type()), 
                       rhs , x , 0.0 , 0.05, 288, write_out );
    out.close();
    dout.close();
}
