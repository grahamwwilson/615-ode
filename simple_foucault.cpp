#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
using namespace std;
using namespace boost::numeric::odeint;
ofstream out, dout;
const double BigOmega = 4.57655e-5;  // In rad/s 
const double wx=M_PI;
const double wy=M_PI;
typedef boost::array< double , 4 > state_type;
// Specify the simple Foucault ODE using equation 4.24 of Blackburn.
// Has only linear terms but allows for an asymmetry in x and y.
// TODO DOUBLE_CHECK the sign convention for BigOmega

// The rhs of x' = f(x,t)
// For [x] = {x, vx, y, vy}
void rhs( const state_type &x , state_type &dxdt , const double t )
{
   dxdt[0] =  x[1];                                 // d/dt(x)  = vx
   dxdt[1] = -wx*wx*x[0] + 2.0*BigOmega*x[3];       // d/dt(vx) = -w^2 x
   dxdt[2] =  x[3];                                 // d/dt(y)  = vy
   dxdt[3] = -wy*wy*x[2] - 2.0*BigOmega*x[1];       // d/dt(vy) = -w^2 y
}

void write_out( const state_type &x , const double t )
{
// should also keep track of energy ...
    out  << t << ' ' << x[0] << ' ' << x[2] 
              << ' ' << x[1] << ' ' << x[3] << endl;
    dout << t << ' ' << x[0] << ' ' << x[2] << endl;
}

typedef runge_kutta_fehlberg78< state_type > stepper_type;

int main()
{
    state_type x = { 1.0 , 0.0, 0.0,  0.0 }; // initial conditions

    out.open("simple_foucault.out");
    out.precision(16);
    dout.open("simple_foucault.dat");
    dout.precision(16);
    out << "Pi: " << M_PI << endl;
    integrate_n_steps( make_controlled(1e-14,1e-14,stepper_type()), 
                       rhs , x , 0.0 , 0.005, 4000, write_out );
    out.close();
    dout.close();
}
