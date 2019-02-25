p.totalmasses = 7;

% linkage masses
p.m1 = 6;    p.m2 = 6;   p.m3 = 6;  p.m4 = 6;  p.m5 = 2;    p.m6 = 2;   p.m7 = 2;
% linkage lengths
p.l1 = 5;    p.l2 = 5;   p.l3 = 5;  p.l4 = 5;  p.l5 = 1;    p.l6 = 1;   p.l7 = 1;
% gravitational constant
p.g  = 10;

%F#'s below are forces in the x direction (e.g. wind)
p.F1 = 0;    p.F2 = 0;   p.F3 = 0;  p.F4 = 0;  p.F5 = 0;    p.F6 = 0;   p.F7 = 0;

theta0    = [  pi/2-.2  ,  pi/2-.1   ,  pi/2+.1  ,  pi/2+.2   ,    pi     ,    pi*.9     ,    pi    ];
thetadot0 = [    0      ,     0      ,    0      ,     0      ,    0     ,    0     ,    0    ];

reltol = 10e-9;
abstol = 10e-9;
tmax   = 30;
timesteps = 10^3;

q0 = [theta0 thetadot0]';

[tarray , qarray] = string_lights_DAEs(p, q0, reltol, abstol, tmax, timesteps);

conservationCheck_string_lights(q0,tarray,qarray,p);

pause(2);
plot_and_animate_lights(p,tarray,qarray);
