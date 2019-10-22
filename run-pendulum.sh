#!/bin/sh
#
# Unix shell script named run-pendulum.sh
#
# Run all 6 ODE integrations.
#
# Parameters are: time-step, out-filename, dat-filename, stepping algorithm
#
# 0 = Euler
# 1 = Second order Runge Kutta (RK2)
# 2 = Fourth order Runge Kutta (RK4)
#

./pendulum 0.01  "euler_0p01.out" "euler_0p01.dat" 0
./pendulum 0.005 "euler_0p005.out" "euler_0p005.dat" 0

./pendulum 0.01  "rk2_0p01.out" "rk2_0p01.dat" 1
./pendulum 0.005 "rk2_0p005.out" "rk2_0p005.dat" 1

./pendulum 0.01  "rk4_0p01.out" "rk4_0p01.dat" 2
./pendulum 0.005 "rk4_0p005.out" "rk4_0p005.dat" 2

exit
