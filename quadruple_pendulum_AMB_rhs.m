function qdot = quadruple_pendulum_AMB_rhs(tspan,q0,flag,p) 
 
    l1 = p.l1; l2 = p.l2; l3 = p.l3; l4 = p.l4; m1 = p.m1; m2 = p.m2; m3 = p.m3; m4 = p.m4; g = p.g; 
    F1 = p.F1; F2 = p.F2; F3 = p.F3; F4 = p.F4; 
    I_G1 = 1/12*m1*l1^2; I_G2 = 1/12*m2*l2^2; I_G3 = 1/12*m3*l3^2; I_G4 = 1/12*m4*l4^2; 
 
    theta1 = q0(1); theta2 = q0(2); theta3 = q0(3); theta4 = q0(4); thetadot1 = q0(5); thetadot2 = q0(6); thetadot3 = q0(7); thetadot4 = q0(8); 
 
    M(1,1) = - I_G1 - (l1^2*m1)/4 - l1^2*m2 - l1^2*m3 - l1^2*m4 - (l1*l2*m2*cos(theta1 - theta2))/2 - l1*l2*m3*cos(theta1 - theta2) - l1*l2*m4*cos(theta1 - theta2) - (l1*l3*m3*cos(theta1 - theta3))/2 - l1*l3*m4*cos(theta1 - theta3) - (l1*l4*m4*cos(theta1 - theta4))/2; 
    M(1,2) = - I_G2 - (l2^2*m2)/4 - l2^2*m3 - l2^2*m4 - (l1*l2*m2*cos(theta1 - theta2))/2 - l1*l2*m3*cos(theta1 - theta2) - l1*l2*m4*cos(theta1 - theta2) - (l2*l3*m3*cos(theta2 - theta3))/2 - l2*l3*m4*cos(theta2 - theta3) - (l2*l4*m4*cos(theta2 - theta4))/2; 
    M(1,3) = - I_G3 - (l3^2*m3)/4 - l3^2*m4 - (l1*l3*m3*cos(theta1 - theta3))/2 - l1*l3*m4*cos(theta1 - theta3) - (l2*l3*m3*cos(theta2 - theta3))/2 - l2*l3*m4*cos(theta2 - theta3) - (l3*l4*m4*cos(theta3 - theta4))/2; 
    M(1,4) = - I_G4 - (l4^2*m4)/4 - (l1*l4*m4*cos(theta1 - theta4))/2 - (l2*l4*m4*cos(theta2 - theta4))/2 - (l3*l4*m4*cos(theta3 - theta4))/2; 
    M(2,1) = -(l1*(l2*m2*cos(theta1 - theta2) + 2*l2*m3*cos(theta1 - theta2) + 2*l2*m4*cos(theta1 - theta2) + l3*m3*cos(theta1 - theta3) + 2*l3*m4*cos(theta1 - theta3) + l4*m4*cos(theta1 - theta4)))/2; 
    M(2,2) = - I_G2 - (l2^2*m2)/4 - l2^2*m3 - l2^2*m4 - (l2*l3*m3*cos(theta2 - theta3))/2 - l2*l3*m4*cos(theta2 - theta3) - (l2*l4*m4*cos(theta2 - theta4))/2; 
    M(2,3) = - I_G3 - (l3^2*m3)/4 - l3^2*m4 - (l2*l3*m3*cos(theta2 - theta3))/2 - l2*l3*m4*cos(theta2 - theta3) - (l3*l4*m4*cos(theta3 - theta4))/2; 
    M(2,4) = - I_G4 - (l4^2*m4)/4 - (l2*l4*m4*cos(theta2 - theta4))/2 - (l3*l4*m4*cos(theta3 - theta4))/2; 
    M(3,1) = -(l1*(l3*m3*cos(theta1 - theta3) + 2*l3*m4*cos(theta1 - theta3) + l4*m4*cos(theta1 - theta4)))/2; 
    M(3,2) = -(l2*(l3*m3*cos(theta2 - theta3) + 2*l3*m4*cos(theta2 - theta3) + l4*m4*cos(theta2 - theta4)))/2; 
    M(3,3) = - I_G3 - (l3^2*m3)/4 - l3^2*m4 - (l3*l4*m4*cos(theta3 - theta4))/2; 
    M(3,4) = - I_G4 - (l4^2*m4)/4 - (l3*l4*m4*cos(theta3 - theta4))/2; 
    M(4,1) = -(l1*l4*m4*cos(theta1 - theta4))/2; 
    M(4,2) = -(l2*l4*m4*cos(theta2 - theta4))/2; 
    M(4,3) = -(l3*l4*m4*cos(theta3 - theta4))/2; 
    M(4,4) = - I_G4 - (l4^2*m4)/4; 
    rhs(1,1) = m3*(l1*sin(theta1) + l2*sin(theta2) + (l3*sin(theta3))/2)*(l1*thetadot1^2*cos(theta1) + l2*thetadot2^2*cos(theta2) + (l3*thetadot3^2*cos(theta3))/2) - F4*(l1*cos(theta1) + l2*cos(theta2) + l3*cos(theta3) + (l4*cos(theta4))/2) - F3*(l1*cos(theta1) + l2*cos(theta2) + (l3*cos(theta3))/2) - m3*(l1*cos(theta1) + l2*cos(theta2) + (l3*cos(theta3))/2)*(l1*thetadot1^2*sin(theta1) + l2*thetadot2^2*sin(theta2) + (l3*thetadot3^2*sin(theta3))/2) - F2*(l1*cos(theta1) + (l2*cos(theta2))/2) + m4*(l1*sin(theta1) + l2*sin(theta2) + l3*sin(theta3) + (l4*sin(theta4))/2)*(l1*thetadot1^2*cos(theta1) + l2*thetadot2^2*cos(theta2) + l3*thetadot3^2*cos(theta3) + (l4*thetadot4^2*cos(theta4))/2) - m4*(l1*thetadot1^2*sin(theta1) + l2*thetadot2^2*sin(theta2) + l3*thetadot3^2*sin(theta3) + (l4*thetadot4^2*sin(theta4))/2)*(l1*cos(theta1) + l2*cos(theta2) + l3*cos(theta3) + (l4*cos(theta4))/2) + g*m4*(l1*sin(theta1) + l2*sin(theta2) + l3*sin(theta3) + (l4*sin(theta4))/2) + g*m3*(l1*sin(theta1) + l2*sin(theta2) + (l3*sin(theta3))/2) + m2*(l1*thetadot1^2*cos(theta1) + (l2*thetadot2^2*cos(theta2))/2)*(l1*sin(theta1) + (l2*sin(theta2))/2) - m2*(l1*thetadot1^2*sin(theta1) + (l2*thetadot2^2*sin(theta2))/2)*(l1*cos(theta1) + (l2*cos(theta2))/2) - (F1*l1*cos(theta1))/2 + g*m2*(l1*sin(theta1) + (l2*sin(theta2))/2) + (g*l1*m1*sin(theta1))/2; 
    rhs(2,1) = m3*(l2*sin(theta2) + (l3*sin(theta3))/2)*(l1*thetadot1^2*cos(theta1) + l2*thetadot2^2*cos(theta2) + (l3*thetadot3^2*cos(theta3))/2) - F4*(l2*cos(theta2) + l3*cos(theta3) + (l4*cos(theta4))/2) - F3*(l2*cos(theta2) + (l3*cos(theta3))/2) - m3*(l2*cos(theta2) + (l3*cos(theta3))/2)*(l1*thetadot1^2*sin(theta1) + l2*thetadot2^2*sin(theta2) + (l3*thetadot3^2*sin(theta3))/2) + m4*(l2*sin(theta2) + l3*sin(theta3) + (l4*sin(theta4))/2)*(l1*thetadot1^2*cos(theta1) + l2*thetadot2^2*cos(theta2) + l3*thetadot3^2*cos(theta3) + (l4*thetadot4^2*cos(theta4))/2) - m4*(l2*cos(theta2) + l3*cos(theta3) + (l4*cos(theta4))/2)*(l1*thetadot1^2*sin(theta1) + l2*thetadot2^2*sin(theta2) + l3*thetadot3^2*sin(theta3) + (l4*thetadot4^2*sin(theta4))/2) + g*m4*(l2*sin(theta2) + l3*sin(theta3) + (l4*sin(theta4))/2) - (F2*l2*cos(theta2))/2 + g*m3*(l2*sin(theta2) + (l3*sin(theta3))/2) + (l2*m2*sin(theta2)*(l1*thetadot1^2*cos(theta1) + (l2*thetadot2^2*cos(theta2))/2))/2 - (l2*m2*cos(theta2)*(l1*thetadot1^2*sin(theta1) + (l2*thetadot2^2*sin(theta2))/2))/2 + (g*l2*m2*sin(theta2))/2; 
    rhs(3,1) = m4*(l3*sin(theta3) + (l4*sin(theta4))/2)*(l1*thetadot1^2*cos(theta1) + l2*thetadot2^2*cos(theta2) + l3*thetadot3^2*cos(theta3) + (l4*thetadot4^2*cos(theta4))/2) - F4*(l3*cos(theta3) + (l4*cos(theta4))/2) - (F3*l3*cos(theta3))/2 - m4*(l3*cos(theta3) + (l4*cos(theta4))/2)*(l1*thetadot1^2*sin(theta1) + l2*thetadot2^2*sin(theta2) + l3*thetadot3^2*sin(theta3) + (l4*thetadot4^2*sin(theta4))/2) + g*m4*(l3*sin(theta3) + (l4*sin(theta4))/2) - (l3*m3*cos(theta3)*(l1*thetadot1^2*sin(theta1) + l2*thetadot2^2*sin(theta2) + (l3*thetadot3^2*sin(theta3))/2))/2 + (g*l3*m3*sin(theta3))/2 + (l3*m3*sin(theta3)*(l1*thetadot1^2*cos(theta1) + l2*thetadot2^2*cos(theta2) + (l3*thetadot3^2*cos(theta3))/2))/2; 
    rhs(4,1) = -(l4*(F4*cos(theta4) - g*m4*sin(theta4) + l1*m4*thetadot1^2*sin(theta1 - theta4) + l2*m4*thetadot2^2*sin(theta2 - theta4) + l3*m4*thetadot3^2*sin(theta3 - theta4)))/2; 
 
    thetadotdot = M\rhs; 
 
    qdot = [thetadot1; thetadot2; thetadot3; thetadot4; thetadotdot ];
 
end 
