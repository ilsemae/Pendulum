This project was completed in 2013 as the final project for the Intermediate Dynamics class at Cornell University, Mechanical Engineering.

This code comprises simulation (in Matlab) of a triple pendulum and a four-bar linkage using three different methods for obtaining the necessary equations. Also included are a quadruple pendulum, a five-bar linkage, and a system that I have named "String of Lights", which is a combination of the five-bar linkage and single pendulums hinged at the unfixed hinges of the five-bar linkage. In all simulations, friction forces are ignored, and the rods are infinitely thin, with the mass distributed evenly throughout their lengths. In some of the simulations, you can add external, constant, nonconservative forces in the horizontal direction.

There are a total of 9 different simulations. All but one can be run using runner.m. Open this file and uncomment one of lines 15, 18, 22, 26, 30, 34, 37, or 41 to run a simulation. To run the last simulation, "string of lights," open and run "runner_lights.m" 

To compare different derivation methods for the triple pendulum or for the "branched pendulum," open and run "lagrange_vs_AMB_vs_DAEs.m." 

Initial conditions for all of these can be set at the top of the files.
