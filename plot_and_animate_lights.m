function plot_and_animate_lights(p,t,q)

    %plot angles vs time
    figure;
    plot(t,q(:,1),'r',t,q(:,2)-4,'b',t,q(:,3)-8,'k');
        if p.totalmasses == 4
            hold on;
            plot(t,q(:,4)-12,'m');
            hold off;
        end
        if p.totalmasses == 7
            hold on;
            plot(t,q(:,4)-12,'m',t,q(:,5)-16,'c',t,q(:,6)-20,'y',t,q(:,7)-24,'g');
            hold off;
        end
    title('rod angles vs time');
    xlabel('time (s)');
    ylabel('angle (radians)');
    
    rod1x = [ zeros(size(q,1),1)    p.l1*sin(q(:,1))];
    rod1y = [ zeros(size(q,1),1)   -p.l1*cos(q(:,1))];
    rod2x = [ rod1x(:,2)            p.l1*sin(q(:,1)) + p.l2*sin(q(:,2))];
    rod2y = [ rod1y(:,2)           -p.l1*cos(q(:,1)) - p.l2*cos(q(:,2))];
    rod3x = [ rod2x(:,2)            p.l1*sin(q(:,1)) + p.l2*sin(q(:,2)) + p.l3*sin(q(:,3))];
    rod3y = [ rod2y(:,2)           -p.l1*cos(q(:,1)) - p.l2*cos(q(:,2)) - p.l3*cos(q(:,3))];
    rod4x = [ rod3x(:,2)            p.l1*sin(q(:,1)) + p.l2*sin(q(:,2)) + p.l3*sin(q(:,3)) + p.l4*sin(q(:,4))];
    rod4y = [ rod3y(:,2)           -p.l1*cos(q(:,1)) - p.l2*cos(q(:,2)) - p.l3*cos(q(:,3)) - p.l4*cos(q(:,4))];
    rod5x = [ rod1x(:,2)            p.l1*sin(q(:,1)) + p.l5*sin(q(:,5))];
    rod5y = [ rod1y(:,2)           -p.l1*cos(q(:,1)) - p.l5*cos(q(:,5))];
    rod6x = [ rod2x(:,2)            p.l1*sin(q(:,1)) + p.l2*sin(q(:,2)) + p.l6*sin(q(:,6))];
    rod6y = [ rod2y(:,2)           -p.l1*cos(q(:,1)) - p.l2*cos(q(:,2)) - p.l6*cos(q(:,6))];
    rod7x = [ rod3x(:,2)            p.l1*sin(q(:,1)) + p.l2*sin(q(:,2)) + p.l3*sin(q(:,3)) + p.l7*sin(q(:,7))];
    rod7y = [ rod3y(:,2)           -p.l1*cos(q(:,1)) - p.l2*cos(q(:,2)) - p.l3*cos(q(:,3)) - p.l7*cos(q(:,7))];
    pause(2);

    %animation
    linkage_width = abs(rod4x(1,2)-rod1x(1,1));
    linkage_height = abs(rod1y(1,1)-rod2y(1,2));
    figure;
    for j = 1:size(t,1)    
        pause(0.01);
        plot(rod1x(j,:),rod1y(j,:),'b',rod2x(j,:),rod2y(j,:),'b',rod3x(j,:),rod3y(j,:),'b',rod4x(j,:),rod4y(j,:),'b',...
             rod5x(j,:),rod5y(j,:),'y',rod6x(j,:),rod6y(j,:),'y',rod7x(j,:),rod7y(j,:),'y','linewidth',4);
        axis equal;
        axis([-5 linkage_width+5 -linkage_height-5 linkage_height+5]);
        if j == 1 pause(2); end
    end

end