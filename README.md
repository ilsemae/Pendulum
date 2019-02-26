There are a total of 9 different simulations. All but one can be run using 
runner.m. Open this file and uncomment one of lines 15, 18, 22, 26, 30, 
34, 37, or 41 to run a simulation. To run the last simulation, "string of 
lights," open and run "runner_lights.m" 

To compare different derivation methods for the triple pendulum or for the 
"branched pendulum," open and run "lagrange_vs_AMB_vs_DAEs.m." 

Initial conditions for all of these can be set at the top of the files.

# Description

## Introduction

This code comprises  simulation (in Matlab) of a triple pendulum and a four-bar linkage using three different methods for obtaining the necessary equations. Also included are a quadruple pendulum, a five-bar linkage, and a system that I have named "String of Lights", which is a combination of the five-bar linkage and single pendulums hinged at the moving hinges of the five-bar linkage. In all simulations, friction forces are ignored, and the rods are infinitely thin, with the mass distributed evenly throughout their lengths. In some of the simulations, you can add external, constant, nonconservative forces in the horizontal direction.

## Triple Pendulum

A triple pendulum is a single pendulum with another pendulum hinged at its end, which in turn has another pendulum hinged at its end. See Appendix 2.01 for a drawing of our triple pendulum setup. In order to find the motion of this system given initial positions and velocities, we need to know the accelerations of at least the angular positions of each linkage. In simpler systems, (such as systems with point masses) we can often use linear momentum balance to find equations of motion. In this case, however, doing so would leave us
with unknown values for the tension forces at the hinges. We can, however, use angular momentum balance, which eliminates these terms.
