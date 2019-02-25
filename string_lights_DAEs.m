% four-bar-linkage using DAEs

function [t , q] = string_lights_DAEs(p,q0,reltol,abstol,tmax,timesteps)
    
    options = odeset('reltol',reltol,'abstol',abstol);
    [t , q] = ode45(@four_bar_linkage_DAEs_rhs,linspace(0,tmax,timesteps),q0,options,p);

end

function qdot = four_bar_linkage_DAEs_rhs(t,q0,p)

    i = [1 0 0]; j = [0 1 0]; k = [0 0 1];

    theta1 = q0(1);
    theta2 = q0(2);
    theta3 = q0(3);
    theta4 = q0(4);
    theta5 = q0(5);
    theta6 = q0(6);
    theta7 = q0(7);
    thetadot1 = q0(8);
    thetadot2 = q0(9);
    thetadot3 = q0(10);  
    thetadot4 = q0(11);  
    thetadot5 = q0(12);
    thetadot6 = q0(13);  
    thetadot7 = q0(14);  
    
    I1 =  1/12*p.m1*p.l1^2;
    I2 =  1/12*p.m2*p.l2^2;
    I3 =  1/12*p.m3*p.l3^2;
    I4 =  1/12*p.m4*p.l4^2;
    I5 =  1/12*p.m5*p.l5^2;
    I6 =  1/12*p.m6*p.l6^2;
    I7 =  1/12*p.m7*p.l7^2;
    
    r_0_to_G1  = p.l1/2 * (sin(theta1)*i - cos(theta1)*j);  r_G1_to_E1 =  r_0_to_G1;
    r_E1_to_G2 = p.l2/2 * (sin(theta2)*i - cos(theta2)*j);  r_G2_to_E2 =  r_E1_to_G2;
    r_E2_to_G3 = p.l3/2 * (sin(theta3)*i - cos(theta3)*j);  r_G3_to_E3 =  r_E2_to_G3;
    r_E3_to_G4 = p.l4/2 * (sin(theta4)*i - cos(theta4)*j);  r_G4_to_E4 =  r_E3_to_G4;
    r_E1_to_G5 = p.l5/2 * (sin(theta5)*i - cos(theta5)*j);  
    r_E2_to_G6 = p.l6/2 * (sin(theta6)*i - cos(theta6)*j);  
    r_E3_to_G7 = p.l7/2 * (sin(theta7)*i - cos(theta7)*j);  
    
    M  = diag([p.m1 , p.m1 , p.m2 , p.m2 , p.m3 , p.m3 , p.m4 , p.m4 , p.m5 , p.m5 , p.m6 , p.m6 , p.m7 , p.m7 , I1 , I2 , I3 , I4 , I5 , I6 , I7]);
    
    %      T0x              T0y            T1x             T1y               T2x            T2y              T3x              T3y              T4x              T4y              T5x            T5y              T6x              T6y              T7x              T7y
    J1 = [  1        ,       0       ,     -1        ,      0         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,      -1       ,      0         ,      0        ,       0         ,      0        ,       0         ;...    LMB of bar1
            0        ,       1       ,      0        ,     -1         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,       0       ,     -1         ,      0        ,       0         ,      0        ,       0         ;...    LMB of bar1
            0        ,       0       ,      1        ,      0         ,      -1       ,      0         ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,     -1        ,       0         ,      0        ,       0         ;...    LMB of bar2
            0        ,       0       ,      0        ,      1         ,       0       ,     -1         ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,      0        ,      -1         ,      0        ,       0         ;...    LMB of bar2
            0        ,       0       ,      0        ,      0         ,       1       ,      0         ,     -1        ,       0         ,      0        ,       0        ,       0       ,      0         ,      0        ,       0         ,     -1        ,       0         ;...    LMB of bar3
            0        ,       0       ,      0        ,      0         ,       0       ,      1         ,      0        ,      -1         ,      0        ,       0        ,       0       ,      0         ,      0        ,       0         ,      0        ,      -1         ;...    LMB of bar3
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      1        ,       0         ,     -1        ,       0        ,       0       ,      0         ,      0        ,       0         ,      0        ,       0         ;...    LMB of bar4
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      0        ,       1         ,      0        ,      -1        ,       0       ,      0         ,      0        ,       0         ,      0        ,       0         ;...    LMB of bar4
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,       1       ,      0         ,      0        ,       0         ,      0        ,       0         ;...    LMB of bar5
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,       0       ,      1         ,      0        ,       0         ,      0        ,       0         ;...    LMB of bar5
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,      1        ,       0         ,      0        ,       0         ;...    LMB of bar6
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,      0        ,       1         ,      0        ,       0         ;...    LMB of bar6
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,      0        ,       0         ,      1        ,       0         ;...    LMB of bar7
            0        ,       0       ,      0        ,      0         ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,      0        ,       0         ,      0        ,       1         ;...    LMB of bar7
        r_0_to_G1(2) , -r_0_to_G1(1) , r_G1_to_E1(2) , -r_G1_to_E1(1) ,       0       ,      0         ,      0        ,       0         ,      0        ,       0        ,  r_G1_to_E1(2), -r_G1_to_E1(1) ,      0        ,       0         ,      0        ,       0         ;...    AMB of bar1
            0        ,       0       , r_E1_to_G2(2) , -r_E1_to_G2(1) , r_G2_to_E2(2) , -r_G2_to_E2(1) ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,  r_G2_to_E2(2), -r_G2_to_E2(1)  ,      0        ,       0         ;...    AMB of bar2    
            0        ,       0       ,      0        ,      0         , r_E2_to_G3(2) , -r_E2_to_G3(1) , r_G3_to_E3(2) , -r_G3_to_E3(1)  ,      0        ,       0        ,       0       ,      0         ,      0        ,       0         ,  r_G3_to_E3(2), -r_G3_to_E3(1)  ;...    AMB of bar3
            0        ,       0       ,      0        ,      0         ,      0        ,       0        , r_E3_to_G4(2) , -r_E3_to_G4(1)  , r_G4_to_E4(2) , -r_G4_to_E4(1) ,       0       ,      0         ,      0        ,       0         ,      0        ,       0         ;...    AMB of bar4
            0        ,       0       ,      0        ,      0         ,      0        ,       0        ,      0        ,       0         ,      0        ,       0        ,  r_E1_to_G5(2), -r_E1_to_G5(1) ,      0        ,       0         ,      0        ,       0         ;...    AMB of bar5
            0        ,       0       ,      0        ,      0         ,      0        ,       0        ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,  r_E2_to_G6(2), -r_E2_to_G6(1)  ,      0        ,       0         ;...    AMB of bar6
            0        ,       0       ,      0        ,      0         ,      0        ,       0        ,      0        ,       0         ,      0        ,       0        ,       0       ,      0         ,      0        ,       0         ,  r_E3_to_G7(2), -r_E3_to_G7(1)  ];...   AMB of bar7
            
    %     xddot1  yddot1   xddot2   yddot2    xddot3   yddot3   xddot4   yddot4    xddot5   yddot5    xddot6   yddot6   xddot7    yddot7     thetaddot1      thetaddot2        thetaddot3      thetaddot4        thetaddot5     thetaddot6      thetaddot7    
    J2 = [  1   ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,  r_0_to_G1(2)  ,      0         ,       0        ,       0        ,      0       ,       0        ,       0       ;... constraint on hinge0
            0   ,   1    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    , -r_0_to_G1(1)  ,      0         ,       0        ,       0        ,      0       ,       0        ,       0       ;... constraint on hinge0
           -1   ,   0    ,   1   ,    0    ,    0    ,   0    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,  r_G1_to_E1(2) ,  r_E1_to_G2(2) ,       0        ,       0        ,      0       ,       0        ,       0       ;... constraint on hinge1
            0   ,  -1    ,   0   ,    1    ,    0    ,   0    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    , -r_G1_to_E1(1) , -r_E1_to_G2(1) ,       0        ,       0        ,      0       ,       0        ,       0       ;... constraint on hinge1
            0   ,   0    ,  -1   ,    0    ,    1    ,   0    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,       0        ,  r_G2_to_E2(2) ,  r_E2_to_G3(2) ,       0        ,      0       ,       0        ,       0       ;... constraint on hinge2 
            0   ,   0    ,   0   ,   -1    ,    0    ,   1    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,       0        , -r_G2_to_E2(1) , -r_E2_to_G3(1) ,       0        ,      0       ,       0        ,       0       ;... constraint on hinge2
            0   ,   0    ,   0   ,    0    ,   -1    ,   0    ,    1    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,       0        ,      0         ,  r_G3_to_E3(2) ,  r_E3_to_G4(2) ,      0       ,       0        ,       0       ;... constraint on hinge3
            0   ,   0    ,   0   ,    0    ,    0    ,  -1    ,    0    ,   1    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,       0        ,      0         , -r_G3_to_E3(1) , -r_E3_to_G4(1) ,      0       ,       0        ,       0       ;... constraint on hinge3
            0   ,   0    ,   0   ,    0    ,    0    ,   0    ,   -1    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,       0        ,      0         ,       0        ,  r_G4_to_E4(2) ,      0       ,       0        ,       0       ;... constraint on hinge4
            0   ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,  -1    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,       0        ,      0         ,       0        , -r_G4_to_E4(1) ,      0       ,       0        ,       0       ;... constraint on hinge4
           -1   ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,   1   ,    0    ,    0    ,   0    ,    0    ,   0    ,  r_G1_to_E1(2) ,      0         ,       0        ,       0        , r_E1_to_G5(2),       0        ,       0       ;... constraint on hinge2b 
            0   ,  -1    ,   0   ,    0    ,    0    ,   0    ,    0    ,   0    ,   0   ,    1    ,    0    ,   0    ,    0    ,   0    , -r_G1_to_E1(1) ,      0         ,       0        ,       0        ,-r_E1_to_G5(1),       0        ,       0       ;... constraint on hinge2b
            0   ,   0    ,  -1   ,    0    ,    0    ,   0    ,    0    ,   0    ,   0   ,    0    ,    1    ,   0    ,    0    ,   0    ,       0        ,  r_G2_to_E2(2) ,       0        ,       0        ,      0       ,  r_E2_to_G6(2) ,       0       ;... constraint on hinge3b
            0   ,   0    ,   0   ,   -1    ,    0    ,   0    ,    0    ,   0    ,   0   ,    0    ,    0    ,   1    ,    0    ,   0    ,       0        , -r_G2_to_E2(1) ,       0        ,       0        ,      0       , -r_E2_to_G6(1) ,       0       ;... constraint on hinge3b
            0   ,   0    ,   0   ,    0    ,   -1    ,   0    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    1    ,   0    ,       0        ,      0         ,  r_G3_to_E3(2) ,       0        ,      0       ,       0        , r_E3_to_G7(2) ;... constraint on hinge4b
            0   ,   0    ,   0   ,    0    ,    0    ,  -1    ,    0    ,   0    ,   0   ,    0    ,    0    ,   0    ,    0    ,   1    ,       0        ,      0         , -r_G3_to_E3(1) ,       0        ,      0       ,       0        ,-r_E3_to_G7(1) ];... constraint on hinge4b

    B =  [  ... LMB
             p.F1                          ;...
            -p.m1*p.g                      ;...
             p.F2                          ;...
            -p.m2*p.g                      ;...
             p.F3                          ;...
            -p.m3*p.g                      ;...
             p.F4                          ;...
            -p.m4*p.g                      ;...
             p.F5                          ;...
            -p.m5*p.g                      ;...
             p.F6                          ;...
            -p.m6*p.g                      ;...
             p.F7                          ;...
            -p.m7*p.g                      ;...
            ...
            ... AMB
             0                             ;...
             0                             ;...
             0                             ;...
             0                             ;...
             0                             ;...
             0                             ;...
             0                             ;...
            ...  
            ... contraints
             -thetadot1^2*r_0_to_G1(1)     ;...
             -thetadot1^2*r_0_to_G1(2)     ;...
             -thetadot1^2*r_G1_to_E1(1) - thetadot2^2*r_E1_to_G2(1)      ;...
             -thetadot1^2*r_G1_to_E1(2) - thetadot2^2*r_E1_to_G2(2)      ;...
             -thetadot2^2*r_G2_to_E2(1) - thetadot3^2*r_E2_to_G3(1)      ;...
             -thetadot2^2*r_G2_to_E2(2) - thetadot3^2*r_E2_to_G3(2)      ;...
             -thetadot3^2*r_G3_to_E3(1) - thetadot4^2*r_E3_to_G4(1)      ;...
             -thetadot3^2*r_G3_to_E3(2) - thetadot4^2*r_E3_to_G4(2)      ;...
             -thetadot4^2*r_G4_to_E4(1)    ;...
             -thetadot4^2*r_G4_to_E4(2)    ;...
             -thetadot1^2*r_G1_to_E1(1) - thetadot5^2*r_E1_to_G5(1)      ;...
             -thetadot1^2*r_G1_to_E1(2) - thetadot5^2*r_E1_to_G5(2)      ;...
             -thetadot2^2*r_G2_to_E2(1) - thetadot6^2*r_E2_to_G6(1)      ;...
             -thetadot2^2*r_G2_to_E2(2) - thetadot6^2*r_E2_to_G6(2)      ;...
             -thetadot3^2*r_G3_to_E3(1) - thetadot7^2*r_E3_to_G7(1)      ;...
             -thetadot3^2*r_G3_to_E3(2) - thetadot7^2*r_E3_to_G7(2)      ;...
             ];
    
    A = [ M  , J1         ;...
          J2 , zeros(16,16) ];

    Z = A\B;
    
    qdot = [thetadot1; thetadot2; thetadot3 ; thetadot4 ; thetadot5; thetadot6 ; thetadot7 ; Z(15:21)];

end