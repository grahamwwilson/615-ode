#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <cmath>
using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::math;
ofstream out, dout;
const double G = 9.81;  // g in m s^-2 
const double L = 1.50;  // l in m
const double THETA0 = M_PI/3.0;
typedef boost::array< double , 2 > state_type;
double fEnergy( const state_type &x );

// The rhs of x' = f(x,t)
// For [x] = {theta, omega}
void rhs( const state_type &x , state_type &dxdt , const double t )
{
    dxdt[0] =  x[1];                // d/dt(theta)  = omega
    dxdt[1] = -(G/L)*sin(x[0]);     // d/dt(omega) = -g/l sin(theta)
}

void write_out( const state_type &x , const double t )
{
// Also keep track of energy
    out  << t << " " << x[0] << " " << x[1] << " " << fEnergy(x) << endl; 
    dout << t << " " << x[0] << " " << x[1] << " " << fEnergy(x) << endl;
}

double fEnergy( const state_type &x )
{
    double m = 1.0;   // set unit mass - not relevant 
    double energy0 = m*G*L*(1.0 - cos(THETA0) );
    double energy = 0.5*m*L*L*x[1]*x[1] + m*G*L*(1.0-cos(x[0]));
    return (-1.0 + energy/energy0);
}

typedef runge_kutta_fehlberg78< state_type > stepper_type;

int main()
{
    state_type x = { THETA0 , 0.0 }; // initial conditions

    double w0 = sqrt(G/L);
    double T0 = 2.0*M_PI/w0;
    double k = sin(THETA0/2.0);
    double T = T0*(2.0/M_PI)*ellint_1(k);

    cout << "THETA0 = " << THETA0 << endl;
    cout << "w0     = " << w0 << endl;
    cout << "T0     = " << T0 << endl;
    cout << "T      = " << T  << endl;
    cout << "T/T0   = " << T/T0 << endl;

    out.open("simple.out");
    out.precision(16);
    dout.open("simple.dat");
    dout.precision(16);
    integrate_n_steps( make_controlled(1e-14,1e-14,stepper_type()), 
                       rhs , x , 0.0 , 0.005, 4000, write_out );
    out.close();
    dout.close();
}
