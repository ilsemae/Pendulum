function qdot = branching_pendulum_lagrange_rhs(tspan,q0,flag,p) 
 
    l1 = p.l1; l2 = p.l2; l3 = p.l3; m1 = p.m1; m2 = p.m2; m3 = p.m3; g = p.g; 
    IG1 = 1/12*m1*l1^2; IG2 = 1/12*m2*l2^2; IG3 = 1/12*m3*l3^2; 
 
    theta1 = q0(1); theta2 = q0(2); theta3 = q0(3); thetadot1 = q0(4); thetadot2 = q0(5); thetadot3 = q0(6); 
 
    M(1,1) = - IG1 - (l1^2*m1)/4 - l1^2*m2 - l1^2*m3; 
    M(1,2) = -(l1*l2*m2*cos(theta1 - theta2))/2; 
    M(1,3) = -(l1*l3*m3*cos(theta1 - theta3))/2; 
    M(2,1) = -(l1*l2*m2*cos(theta1 - theta2))/2; 
    M(2,2) = - IG2 - (l2^2*m2)/4; 
    M(2,3) = 0; 
    M(3,1) = -(l1*l3*m3*cos(theta1 - theta3))/2; 
    M(3,2) = 0; 
    M(3,3) = - IG3 - (l3^2*m3)/4; 
    rhs(1,1) = (l1*(g*m1*sin(theta1) + 2*g*m2*sin(theta1) + 2*g*m3*sin(theta1) + l2*m2*thetadot2^2*sin(theta1 - theta2) + l3*m3*thetadot3^2*sin(theta1 - theta3)))/2; 
    rhs(2,1) = (l2*m2*(g*sin(theta2) - l1*thetadot1^2*sin(theta1 - theta2)))/2; 
    rhs(3,1) = (l3*m3*(g*sin(theta3) - l1*thetadot1^2*sin(theta1 - theta3)))/2; 
 
    thetadotdot = M\rhs; 
    qdot = [thetadot1; thetadot2; thetadot3; thetadotdot ];
 
end 
