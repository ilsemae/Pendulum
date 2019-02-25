function EOMs_derive_triple_pendulum_lagrange() 

    syms theta1 theta2 theta3 thetadot1 thetadot2 thetadot3 thetadotdot1 thetadotdot2 thetadotdot3 real;
    syms IG1 IG2 IG3 m1 m2 m3 l1 l2 l3 g real;
    
    i = [1 0 0]; j = [0 1 0]; k = [0 0 1];
    
    r_E1_0  = l1  * (sin(theta1)*i - cos(theta1)*j);
    r_G1_0  = r_E1_0/2;
    r_E2_E1 = l2  * (sin(theta2)*i - cos(theta2)*j);
    r_G2_E1 = r_E2_E1/2;
    r_E3_E2 = l3  * (sin(theta3)*i - cos(theta3)*j);
    r_G3_E2 = r_E3_E2/2;
     
    v_G1_0 = cross(thetadot1*k, r_G1_0);
    v_G2_0 = cross(thetadot1*k, r_E1_0)  +  cross(thetadot2*k, r_G2_E1);
    v_G3_0 = cross(thetadot1*k, r_E1_0)  +  cross(thetadot2*k, r_E2_E1)  +  cross(thetadot3*k, r_G3_E2);
    
    KE_1   = 1/2*m1*(v_G1_0*v_G1_0') + 1/2*IG1*thetadot1^2;
    KE_2   = 1/2*m2*(v_G2_0*v_G2_0') + 1/2*IG2*thetadot2^2;
    KE_3   = 1/2*m3*(v_G3_0*v_G3_0') + 1/2*IG3*thetadot3^2;
    KE_tot = KE_1 + KE_2 + KE_3;
    
    PE_1   = m1*g* l1/2*(1-cos(theta1));
    PE_2   = m2*g*(l1*(1-cos(theta1))+l2/2*(1-cos(theta2)));
    PE_3   = m3*g*(l1*(1-cos(theta1))+l2*(1-cos(theta2))+l3/2*(1-cos(theta3)));
    PE_tot = PE_1 + PE_2 + PE_3;
    
    L = KE_tot - PE_tot;
    
    EoM = jacobian(L, [theta1 theta2 theta3])' -  jacobian(jacobian(L,[thetadot1 thetadot2 thetadot3])',[theta1 theta2 theta3 thetadot1 thetadot2 thetadot3])*[thetadot1 thetadot2 thetadot3 thetadotdot1 thetadotdot2 thetadotdot3]';
    M   =  simplify(jacobian(EoM, [thetadotdot1 thetadotdot2 thetadotdot3]));     % mass 'matrix'
    rhs = -simplify(EoM - M*[thetadotdot1 thetadotdot2 thetadotdot3]');  % for use in ODE45 rhs function
    
    rhs_file = fopen('triple_pendulum_lagrange_rhs.m','w');
    fprintf(rhs_file,'function qdot = triple_pendulum_lagrange_rhs(tspan,q0,flag,p) \n');
    fprintf(rhs_file,' \n');
    %unpacking
    fprintf(rhs_file,'    l1 = p.l1; l2 = p.l2; l3 = p.l3; m1 = p.m1; m2 = p.m2; m3 = p.m3; g = p.g; \n');
    fprintf(rhs_file,'    IG1 = 1/12*m1*l1^2; IG2 = 1/12*m2*l2^2; IG3 = 1/12*m3*l3^2; \n');
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
    fprintf(rhs_file,'    thetadotdot = M\\rhs;');
    fprintf(rhs_file,' \n');
    fprintf(rhs_file,'    qdot = [thetadot1; thetadot2; thetadot3; thetadotdot ];\n');
    fprintf(rhs_file,' \n');
    fprintf(rhs_file,'end \n');

end