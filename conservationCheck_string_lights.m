function conservationCheck_string_lights(q0,t,q,p) 

    m = zeros(1,7);
    l = zeros(7,1);
    I = zeros(1,7);
    m(1) = p.m1; m(2) = p.m2; m(3) = p.m3; m(4) = p.m4; m(5) = p.m5; m(6) = p.m6; m(7) = p.m7;
    l(1) = p.l1; l(2) = p.l2; l(3) = p.l3; l(4) = p.l4; l(5) = p.l5; l(6) = p.l6; l(7) = p.l7;
    I(1) = 1/12*m(1)*l(1)^2; I(2) = 1/12*m(2)*l(2)^2; I(3) = 1/12*m(3)*l(3)^2; I(4) = 1/12*m(4)*l(4)^2; I(5) = 1/12*m(5)*l(5)^2; I(6) = 1/12*m(6)*l(6)^2; I(7) = 1/12*m(7)*l(7)^2;

    theta0    =  q0(1:7);
    thetadot0 =  q0(8:14);
    h0        =  [l(1)/2*(1-cos(theta0(1)))   ,...
                  l(1)*(1-cos(theta0(1)))+l(2)/2*(1-cos(theta0(2))) ,...
                  l(1)*(1-cos(theta0(1)))+l(2)*(1-cos(theta0(2)))+l(3)/2*(1-cos(theta0(3))) ,...
                  l(1)*(1-cos(theta0(1)))+l(2)*(1-cos(theta0(2)))+l(3)*(1-cos(theta0(3)))+l(4)/2*(1-cos(theta0(4))) ,...
                  l(1)*(1-cos(theta0(1)))+l(5)/2*(1-cos(theta0(5))) ,...
                  l(1)*(1-cos(theta0(1)))+l(2)*(1-cos(theta0(2)))+l(6)/2*(1-cos(theta0(6))) ,...
                  l(1)*(1-cos(theta0(1)))+l(2)*(1-cos(theta0(2)))+l(3)*(1-cos(theta0(3)))+l(7)/2*(1-cos(theta0(7)))];
    v01       =  l(1)/2*thetadot0(1)* [cos(theta0(1)) , -sin(theta0(1))];
    v02       =  2*v01 + l(2)/2*thetadot0(2)* [cos(theta0(2)) , -sin(theta0(2))];
    v03       =  2*v01 + l(2)  *thetadot0(2)* [cos(theta0(2)) , -sin(theta0(2))] + l(3)/2*thetadot0(3)* [cos(theta0(3)) , -sin(theta0(3))];
    v04       =  2*v01 + l(2)  *thetadot0(2)* [cos(theta0(2)) , -sin(theta0(2))] + l(3)  *thetadot0(3)* [cos(theta0(3)) , -sin(theta0(3))] + l(4)/2*thetadot0(4)* [cos(theta0(4)) , -sin(theta0(4))];
    v05       =  2*v01 + l(5)/2*thetadot0(5)* [cos(theta0(5)) , -sin(theta0(5))];
    v06       =  2*v01 + l(2)  *thetadot0(2)* [cos(theta0(2)) , -sin(theta0(2))] + l(6)/2*thetadot0(6)* [cos(theta0(6)) , -sin(theta0(6))];
    v07       =  2*v01 + l(2)  *thetadot0(2)* [cos(theta0(2)) , -sin(theta0(2))] + l(3)  *thetadot0(3)* [cos(theta0(3)) , -sin(theta0(3))] + l(7)/2*thetadot0(7)* [cos(theta0(7)) , -sin(theta0(7))];

    theta     =  q(:,1:7);
    thetadot  =  q(:,8:14);
    
    h         =  [l(1)/2*(1-cos(theta(:,1)))   ,...
                  l(1)*(1-cos(theta(:,1)))+l(2)/2*(1-cos(theta(:,2)))   ,...
                  l(1)*(1-cos(theta(:,1)))+l(2)  *(1-cos(theta(:,2)))+l(3)/2*(1-cos(theta(:,3)))   ,...
                  l(1)*(1-cos(theta(:,1)))+l(2)  *(1-cos(theta(:,2)))+l(3)  *(1-cos(theta(:,3)))+l(4)/2*(1-cos(theta(:,4)))   ,...
                  l(1)*(1-cos(theta(:,1)))+l(5)/2*(1-cos(theta(:,5)))   ,...
                  l(1)*(1-cos(theta(:,1)))+l(2)  *(1-cos(theta(:,2)))+l(6)/2*(1-cos(theta(:,6)))   ,...
                  l(1)*(1-cos(theta(:,1)))+l(2)  *(1-cos(theta(:,2)))+l(3)  *(1-cos(theta(:,3)))+l(7)/2*(1-cos(theta(:,7)))];
              
    v1 = zeros(size(q,1),2);
    v2 = zeros(size(q,1),2);
    v3 = zeros(size(q,1),2);
    v4 = zeros(size(q,1),2);
    v5 = zeros(size(q,1),2);
    v6 = zeros(size(q,1),2);
    v7 = zeros(size(q,1),2);

    for j=1:size(q,1)
        v1(j,:)   =  l(1)/2*thetadot(j,1)* [cos(theta(j,1)) , -sin(theta(j,1))];
        v2(j,:)   =  2*v1(j,:) + l(2)/2*thetadot(j,2)* [cos(theta(j,2)) , -sin(theta(j,2))];
        v3(j,:)   =  2*v1(j,:) + l(2)  *thetadot(j,2)* [cos(theta(j,2)) , -sin(theta(j,2))] + l(3)/2*thetadot(j,3)* [cos(theta(j,3)) , -sin(theta(j,3))];
        v4(j,:)   =  2*v1(j,:) + l(2)  *thetadot(j,2)* [cos(theta(j,2)) , -sin(theta(j,2))] + l(3)  *thetadot(j,3)* [cos(theta(j,3)) , -sin(theta(j,3))] + l(4)/2*thetadot(j,4)* [cos(theta(j,4)) , -sin(theta(j,4))];
        v5(j,:)   =  2*v1(j,:) + l(5)/2*thetadot(j,5)* [cos(theta(j,5)) , -sin(theta(j,5))];
        v6(j,:)   =  2*v1(j,:) + l(2)  *thetadot(j,2)* [cos(theta(j,2)) , -sin(theta(j,2))] + l(6)/2*thetadot(j,6)* [cos(theta(j,6)) , -sin(theta(j,6))];
        v7(j,:)   =  2*v1(j,:) + l(2)  *thetadot(j,2)* [cos(theta(j,2)) , -sin(theta(j,2))] + l(3)  *thetadot(j,3)* [cos(theta(j,3)) , -sin(theta(j,3))] + l(7)/2*thetadot(j,7)* [cos(theta(j,7)) , -sin(theta(j,7))];
    end

    ke = zeros(size(q,1),1);
    pe = zeros(size(q,1),1);

    ke0 = 1/2*(m(1)*(v01*v01')+I(1)*thetadot0(1)^2 +...COM1
               m(2)*(v02*v02')+I(2)*thetadot0(2)^2 +...COM2
               m(3)*(v03*v03')+I(3)*thetadot0(3)^2 +...COM3
               m(4)*(v04*v04')+I(4)*thetadot0(4)^2 +...COM4
               m(5)*(v05*v05')+I(5)*thetadot0(5)^2 +...COM5
               m(6)*(v06*v06')+I(6)*thetadot0(6)^2 +...COM6
               m(7)*(v07*v07')+I(7)*thetadot0(7)^2 );...COM7
    pe0 = p.g*(m(1)*h0(1) +... COM1
               m(2)*h0(2) +... COM2
               m(3)*h0(3) +... COM3
               m(4)*h0(4) +... COM4
               m(5)*h0(5) +... COM5
               m(6)*h0(6) +... COM6
               m(7)*h0(7) ); % COM7
    for j=1:size(q,1)
        ke(j) = 1/2*(m(1)*(v1(j,:)*v1(j,:)')+I(1)*thetadot(j,1)^2 +... COM1
                     m(2)*(v2(j,:)*v2(j,:)')+I(2)*thetadot(j,2)^2 +... COM2
                	 m(3)*(v3(j,:)*v3(j,:)')+I(3)*thetadot(j,3)^2 +... COM3
                	 m(4)*(v4(j,:)*v4(j,:)')+I(4)*thetadot(j,4)^2 +... COM4
                	 m(5)*(v5(j,:)*v5(j,:)')+I(5)*thetadot(j,5)^2 +... COM5
                	 m(6)*(v6(j,:)*v6(j,:)')+I(6)*thetadot(j,6)^2 +... COM6
                	 m(7)*(v7(j,:)*v7(j,:)')+I(7)*thetadot(j,7)^2 ); % COM7
        pe(j) = p.g*(m(1)*h(j,1) +...   COM1
                     m(2)*h(j,2) +...   COM2
                     m(3)*h(j,3) +...   COM3
                     m(4)*h(j,4) +...   COM4
                     m(5)*h(j,5) +...   COM5
                     m(6)*h(j,6) +...   COM6
                     m(7)*h(j,7) );     % COM7
    end
    
    error_energy = ke+pe-ke0*ones(size(q,1),1)-pe0*ones(size(q,1),1);
    figure;
    plot(t,error_energy)
    xlabel('time (s)');
    ylabel('error in energy');
    title('error in energy vs time');
    %figure;
    %plot(t,pe,'b',t,pe0*ones(size(q,1),1),'b -.')
    %figure;
    %plot(t,ke,'r',t,ke0*ones(size(q,1),1),'r -.')

end