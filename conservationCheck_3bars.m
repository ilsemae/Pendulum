function conservationCheck_3bars(q0,t,q,p) 

    m = zeros(1,3);
    l = zeros(3,1);
    I = zeros(1,3);
    m(1) = p.m1; m(2) = p.m2; m(3) = p.m3;
    l(1) = p.l1; l(2) = p.l2; l(3) = p.l3;
    I(1) = 1/12*m(1)*l(1)^2; I(2) = 1/12*m(2)*l(2)^2; I(3) = 1/12*m(3)*l(3)^2;

    theta0    =  q0(1:3);
    thetadot0 =  q0(4:6);
    h0        =  [l(1)/2*(1-cos(theta0(1)))      l(1)*(1-cos(theta0(1)))+l(2)/2*(1-cos(theta0(2))) l(1)*(1-cos(theta0(1)))+l(2)*(1-cos(theta0(2)))+l(3)/2*(1-cos(theta0(3)))];
    v01       =  l(1)/2*thetadot0(1)* [cos(theta0(1)) , -sin(theta0(1))];
    v02       =  2*v01 + l(2)/2*thetadot0(2)* [cos(theta0(2)) , -sin(theta0(2))];
    v03       =  2*v01 + l(2)*thetadot0(2)* [cos(theta0(2)) , -sin(theta0(2))] + l(3)/2*thetadot0(3)* [cos(theta0(3)) , -sin(theta0(3))];

    theta     =  q(:,1:3);
    thetadot  =  q(:,4:6);
    h         =  [l(1)/2*(1-cos(theta(:,1)))   ,   l(1)*(1-cos(theta(:,1)))+l(2)/2*(1-cos(theta(:,2)))   ,   l(1)*(1-cos(theta(:,1)))+l(2)*(1-cos(theta(:,2)))+l(3)/2*(1-cos(theta(:,3)))];

    v1 = zeros(size(q,1),2);
    v2 = zeros(size(q,1),2);
    v3 = zeros(size(q,1),2);

    for j=1:size(q,1)
        v1(j,:)   =  l(1)/2*thetadot(j,1)* [cos(theta(j,1)) , -sin(theta(j,1))];
        v2(j,:)   =  2*v1(j,:) + l(2)/2*thetadot(j,2)* [cos(theta(j,2)) , -sin(theta(j,2))];
        v3(j,:)   =  2*v1(j,:) + l(2)*thetadot(j,2)* [cos(theta(j,2)) , -sin(theta(j,2))] + l(3)/2*thetadot(j,3)* [cos(theta(j,3)) , -sin(theta(j,3))];
    end

    ke = zeros(size(q,1),1);
    pe = zeros(size(q,1),1);

    ke0 = 1/2*m(1)*(v01*v01')+1/2*I(1)*thetadot0(1)^2 +...COM1
          1/2*m(2)*(v02*v02')+1/2*I(2)*thetadot0(2)^2 +...COM2
          1/2*m(3)*(v03*v03')+1/2*I(3)*thetadot0(3)^2;  % COM3
    pe0 = m(1)*p.g*h0(1) +... COM1
          m(2)*p.g*h0(2) +... COM2
          m(3)*p.g*h0(3) ;  % COM3

    for j=1:size(q,1)
        ke(j) = 1/2*m(1)*(v1(j,:)*v1(j,:)')+1/2*I(1)*thetadot(j,1)^2 +... COM1
                1/2*m(2)*(v2(j,:)*v2(j,:)')+1/2*I(2)*thetadot(j,2)^2 +... COM2
                1/2*m(3)*(v3(j,:)*v3(j,:)')+1/2*I(3)*thetadot(j,3)^2;   % COM3
        pe(j) = m(1)*p.g*h(j,1) +...   COM1
                m(2)*p.g*h(j,2) +...   COM2
                m(3)*p.g*h(j,3);     % COM3
    end
    
    error_energy = ke+pe-ke0*ones(size(q,1),1)-pe0*ones(size(q,1),1);
        
    figure;
    plot(t,error_energy);
    xlabel('time (s)');
    ylabel('error in energy');
    title('error in energy vs time');

end