
% linkage masses
p.m1 = 10;       p.m2 = 10;      p.m3 = 10;       p.m4 = 10;
% linkage lengths
p.l1 = 5;       p.l2 = 5;       p.l3 = 5;        p.l4 = 5;
% gravitational constant
p.g  = 10;
% initial angle and angular speeds
theta0    = [  pi/2  ,   0    ,   pi/4    ,   5*pi/4   ];
thetadot0 = [    0   ,   0    ,   0    ,   0   ];

reltol = 10e-9;
abstol = 10e-9;
tmax   = 50;
timesteps = 10^3;

%triple pendulum via angular momentum balance
%**F1-F3 are external forces in the x direction
%   p.F1 = 0; p.F2 = 0; p.F3 = 0; p.totalmasses = 3; q0 = [theta0(1:3) thetadot0(1:3)]';  [tarray , qarray] =         triple_pendulum_AMB(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_3bars(q0,tarray,qarray,p);      plot_and_animate(p,tarray,qarray);

%triple pendulum via lagrange equations
%    p.totalmasses = 3; q0 = [theta0(1:3) thetadot0(1:3)]';  [tarray , qarray] =    triple_pendulum_lagrange(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_3bars(q0,tarray,qarray,p);   disp('hi');   plot_and_animate(p,tarray,qarray); disp('bye');

%triple pendulum via differential algebraic equations
%**F1-F3 are external forces in the x direction
%   p.F1 = 0; p.F2 = 0; p.F3 = 0; p.totalmasses = 3; q0 = [theta0(1:3) thetadot0(1:3)]';  [tarray , qarray] =        triple_pendulum_DAEs(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_3bars(q0,tarray,qarray,p);      plot_and_animate(p,tarray,qarray);

%quadruple pendulum via AMB
%**F1-F4 are external forces in the x direction
  p.F1 = 0; p.F2 = 0; p.F3 = 0; p.F4 = 0;   p.totalmasses = 4; q0 = [theta0 thetadot0]';            [tarray , qarray] =      quadruple_pendulum_AMB(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_4bars(q0,tarray,qarray,p);      plot_and_animate(p,tarray,qarray);

%4-bar linkage via DAEs
%**F1-F3 are external forces in the x direction
%   p.F1 = 0; p.F2 = 0; p.F3 = 0;   p.totalmasses = 3; q0 = [theta0(1:3) thetadot0(1:3)]';  [tarray , qarray] =       four_bar_linkage_DAEs(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_3bars(q0,tarray,qarray,p);      plot_and_animate(p,tarray,qarray);

%5-bar linkage via DAEs
%**F1-F4 are external forces in the x direction
%   p.F1 = 0; p.F2 = 0; p.F3 = 0; p.F4 = 0;   p.totalmasses = 4; q0 = [theta0 thetadot0]';            [tarray , qarray] =       five_bar_linkage_DAEs(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_4bars(q0,tarray,qarray,p);      plot_and_animate(p,tarray,qarray);

%branching double pendulum via lagrange equations
%   p.totalmasses = 3; q0 = [theta0(1:3) thetadot0(1:3)]';  [tarray , qarray] = branching_pendulum_lagrange(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_branching(q0,tarray,qarray,p);  plot_and_animate_branching(p,tarray,qarray);

%branching double pendulum via DAEs
%**F1-F3 are external forces in the x direction
%   p.F1 = 0; p.F2 = 0; p.F3 = 0;   p.totalmasses = 3; q0 = [theta0(1:3) thetadot0(1:3)]';  [tarray , qarray] =     branching_pendulum_DAEs(p, q0, reltol, abstol, tmax, timesteps); conservationCheck_branching(q0,tarray,qarray,p);  plot_and_animate_branching(p,tarray,qarray);

