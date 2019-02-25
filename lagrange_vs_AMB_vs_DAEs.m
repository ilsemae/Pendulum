
p.m1 = 6;        p.m2 = 6;       p.m3 = 6;        p.m4 = 6;
p.l1 = 4;        p.l2 = 4;       p.l3 = 4;        p.l4 = 5;
%do not change these - lagrange equations do not take these forces.
p.F1 = 0;        p.F2 = 0;       p.F3 = 0;        p.F4 = 0;
p.g  = 10;
theta0    = [    .5   ,    pi/2      ,      pi+.5      ,      0    ];
thetadot0 = [    0    ,    0      ,      0      ,      0    ];

reltol = 10e-9;
abstol = 10e-9;
tmax   = 100;
timesteps = 10^3;

q0 = [theta0(1:3) thetadot0(1:3)]';

%triple pendulum
[t , q_triple_pendulum_AMB] =      triple_pendulum_AMB(p, q0, reltol, abstol, tmax, timesteps);
[t , q_triple_pendulum_DAEs] =     triple_pendulum_DAEs(p, q0, reltol, abstol, tmax, timesteps);
[t , q_triple_pendulum_lagrange] = triple_pendulum_lagrange(p, q0, reltol, abstol, tmax, timesteps);

figure; plot(t,q_triple_pendulum_AMB(:,1) -q_triple_pendulum_DAEs(:,1),    t,q_triple_pendulum_AMB(:,2) -q_triple_pendulum_DAEs(:,2)-5,    t,q_triple_pendulum_AMB(:,3) -q_triple_pendulum_DAEs(:,3)-10);
title('difference between AMB and DAE solutions for triple pendulum');
xlabel('time');
legend('mass 1','mass 2','mass 3');
figure; plot(t,q_triple_pendulum_AMB(:,1) -q_triple_pendulum_lagrange(:,1),t,q_triple_pendulum_AMB(:,2) -q_triple_pendulum_lagrange(:,2)-5,t,q_triple_pendulum_AMB(:,3) -q_triple_pendulum_lagrange(:,3)-10);
title('difference between AMB and Lagrange solutions for triple pendulum');
xlabel('time');
legend('mass 1','mass 2','mass 3');
figure; plot(t,q_triple_pendulum_DAEs(:,1)-q_triple_pendulum_lagrange(:,1),t,q_triple_pendulum_DAEs(:,2)-q_triple_pendulum_lagrange(:,2)-5,t,q_triple_pendulum_DAEs(:,3)-q_triple_pendulum_lagrange(:,3)-10);
title('difference between DAE and Lagrange solutions for triple pendulum');
xlabel('time');
legend('mass 1','mass 2','mass 3');

%branched pendulum
[t , q_branching_pendulum_lagrange] = branching_pendulum_lagrange(p, q0, reltol, abstol, tmax, timesteps);
[t , q_branching_pendulum_DAEs] =     branching_pendulum_DAEs(p, q0, reltol, abstol, tmax, timesteps);

figure; plot(t,q_branching_pendulum_lagrange(:,1)-q_branching_pendulum_DAEs(:,1),t,q_branching_pendulum_lagrange(:,2)-q_branching_pendulum_DAEs(:,2)-5,t,q_branching_pendulum_lagrange(:,3)-q_branching_pendulum_DAEs(:,3)-10);
title('difference between DAE and Lagrange solutions for branched pendulum');
xlabel('time');
legend('mass 1','mass 2','mass 3');
