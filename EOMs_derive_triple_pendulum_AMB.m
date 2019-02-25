function EOMs_derive_triple_pendulum_AMB()

    syms theta1 theta2 theta3 thetadot1 thetadot2 thetadot3 thetadotdot1 thetadotdot2 thetadotdot3 real;
    syms l1 l2 l3 m1 m2 m3 g I_G1 I_G2 I_G3 F1 F2 F3 real;

    i = [1 0 0]; j = [0 1 0]; k = [0 0 1];
    x1 = sin(theta1)*i - cos(theta1)*j; y1 = cos(theta1)*i + sin(theta1)*j;
    x2 = sin(theta2)*i - cos(theta2)*j; y2 = cos(theta2)*i + sin(theta2)*j;
    x3 = sin(theta3)*i - cos(theta3)*j; y3 = cos(theta3)*i + sin(theta3)*j;
    
    r_0_to_E1  = l1 * x1;
    r_E1_to_E2 = l2 * x2;
    r_E2_to_G3 = l3 * x3  / 2;
    r_0_to_G1  = r_0_to_E1 / 2 ;
    r_0_to_G2  = r_0_to_E1 + r_E1_to_E2 / 2 ;
    r_0_to_G3  = r_0_to_E1 + r_E1_to_E2 + r_E2_to_G3 ;
    r_E1_to_G2 = r_E1_to_E2 / 2;
    r_E1_to_G3 = r_E1_to_E2 + r_E2_to_G3;
    
    a_0_to_E1  = l1 * (thetadotdot1 * y1 - thetadot1^2 * x1);
    a_E1_to_E2 = l2 * (thetadotdot2 * y2 - thetadot2^2 * x2);
    a_E2_to_P  = l3 * (thetadotdot3 * y3 - thetadot3^2 * x3);
    a_0_to_G1  = a_0_to_E1 / 2 ;
    a_0_to_G2  = a_0_to_E1 + a_E1_to_E2 / 2 ;
    a_0_to_G3  = a_0_to_E1 + a_E1_to_E2 + a_E2_to_P / 2 ;
    
    sumM_0  = cross(r_0_to_G1, -m1*g*j) + cross(r_0_to_G2, -m2*g*j) + cross(r_0_to_G3, -m3*g*j) + cross(r_0_to_G1, F1*i) + cross(r_0_to_G2,F2*i) + cross(r_0_to_G3,F3*i);
    sumM_E1 = cross(r_E1_to_G2,-m2*g*j) + cross(r_E1_to_G3,-m3*g*j) + cross(r_E1_to_G2,F2*i)    + cross(r_E1_to_G3,F3*i);
    sumM_E2 = cross(r_E2_to_G3,-m3*g*j) + cross(r_E2_to_G3,F3*i);
    Hdot_0  = cross(r_0_to_G1, m1*a_0_to_G1) + I_G1*thetadotdot1*k + cross(r_0_to_G2, m2*a_0_to_G2) + I_G2*thetadotdot2*k + cross(r_0_to_G3,m3*a_0_to_G3) + I_G3*thetadotdot3*k;
    Hdot_E1 = cross(r_E1_to_G2,m2*a_0_to_G2) + I_G2*thetadotdot2*k + cross(r_E1_to_G3,m3*a_0_to_G3) + I_G3*thetadotdot3*k;
    Hdot_E2 = cross(r_E2_to_G3,m3*a_0_to_G3) + I_G3*thetadotdot3*k;

    EOM_0  = sumM_0  - Hdot_0;
    EOM_E1 = sumM_E1 - Hdot_E1;
    EOM_E2 = sumM_E2 - Hdot_E2;
    
    EOM_0_scalar  = dot(EOM_0, k);
    EOM_E1_scalar = dot(EOM_E1,k);
    EOM_E2_scalar = dot(EOM_E2,k);
    
    rhs_0                       = subs(EOM_0_scalar, [thetadotdot1 thetadotdot2 thetadotdot3],[0 0 0]);
    thetadotdot1_coefficient_0  = subs(EOM_0_scalar, [thetadotdot1 thetadotdot2 thetadotdot3],[1 0 0])-rhs_0;
    thetadotdot2_coefficient_0  = subs(EOM_0_scalar, [thetadotdot1 thetadotdot2 thetadotdot3],[0 1 0])-rhs_0;
    thetadotdot3_coefficient_0  = subs(EOM_0_scalar, [thetadotdot1 thetadotdot2 thetadotdot3],[0 0 1])-rhs_0;
    rhs_E1                      = subs(EOM_E1_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[0 0 0]);
    thetadotdot1_coefficient_E1 = subs(EOM_E1_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[1 0 0])-rhs_E1;
    thetadotdot2_coefficient_E1 = subs(EOM_E1_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[0 1 0])-rhs_E1;
    thetadotdot3_coefficient_E1 = subs(EOM_E1_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[0 0 1])-rhs_E1;
    rhs_E2                      = subs(EOM_E2_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[0 0 0]);
    thetadotdot1_coefficient_E2 = subs(EOM_E2_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[1 0 0])-rhs_E2;
    thetadotdot2_coefficient_E2 = subs(EOM_E2_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[0 1 0])-rhs_E2;
    thetadotdot3_coefficient_E2 = subs(EOM_E2_scalar,[thetadotdot1 thetadotdot2 thetadotdot3],[0 0 1])-rhs_E2;
    
    rhs = simplify(-[rhs_0 ; rhs_E1 ; rhs_E2]); %negative because we are pulling it to other side of equation
    M   = simplify( [thetadotdot1_coefficient_0  thetadotdot2_coefficient_0  thetadotdot3_coefficient_0 ;...
                     thetadotdot1_coefficient_E1 thetadotdot2_coefficient_E1 thetadotdot3_coefficient_E1 ;...
                     thetadotdot1_coefficient_E2 thetadotdot2_coefficient_E2 thetadotdot3_coefficient_E2 ]);

    rhs_file = fopen('triple_pendulum_AMB_rhs.m','w');
    fprintf(rhs_file,'function qdot = triple_pendulum_AMB_rhs(tspan,q0,flag,p) \n');
    fprintf(rhs_file,' \n');
    %unpacking
    fprintf(rhs_file,'    l1 = p.l1; l2 = p.l2; l3 = p.l3; m1 = p.m1; m2 = p.m2; m3 = p.m3; g = p.g; \n');
    fprintf(rhs_file,'    F1 = p.F1; F2 = p.F2; F3 = p.F3; \n');
    fprintf(rhs_file,'    I_G1 = 1/12*m1*l1^2; I_G2 = 1/12*m2*l2^2; I_G3 = 1/12*m3*l3^2; \n');
    fprintf(rhs_file,' \n');
    fprintf(rhs_file,'    theta1 = q0(1); theta2 = q0(2); theta3 = q0(3); thetadot1 = q0(4); thetadot2 = q0(5); thetadot3 = q0(6); \n');
    fprintf(rhs_file,' \n');
    for i = 1:size(M,1)
        for j = 1:size(M,2)
            fprintf(rhs_file,['    M(' num2str(i) ',' num2str(j) ') = ' char(M(i,j)) '; \n']);
        end
    end
    for i = 1:size(rhs,1)
        for j = 1:size(rhs,2)
            fprintf(rhs_file,['    rhs(' num2str(i) ',' num2str(j) ') = ' char(rhs(i,j)) '; \n']);
        end
    end
    fprintf(rhs_file,' \n');
    fprintf(rhs_file,'    thetadotdot = M\\rhs; \n');
    fprintf(rhs_file,' \n');
    fprintf(rhs_file,'    qdot = [thetadot1; thetadot2; thetadot3; thetadotdot ];\n');
    fprintf(rhs_file,' \n');
    fprintf(rhs_file,'end \n');

end