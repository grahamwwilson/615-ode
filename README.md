# 615-ode
Investigate ODE methods

Currently there are 3 problems implemented

0. pendulum.cpp
Simple finite angle pendulum using explicit Euler, RK2, RK4

1. simple.cpp 
The same simple finite angle pendulum problem implemented 
using boost ODE methods

2. simple_foucault.cpp
Foucault pendulum using boost ODE methods

3. golf.cpp
Golf ball trajectory with drag and spin using boost ODE methods

There is also a reader.py file that can be used to make ROOT TTrees 
from the output *.dat file of the first program (pendulum.cpp).
