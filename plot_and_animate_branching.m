function plot_and_animate_branching(p,t,q)

    %plot angles vs time
    figure;
    plot(t,q(:,1),'r',t,q(:,2)-4,'b',t,q(:,3)-8,'g');
        if p.totalmasses == 4
            hold on;
            plot(t,q(:,4)-12,'m');
            hold off;
        end
    title('rod angles vs time');
    xlabel('time (s)');
    ylabel('angle (radians)');
    
    rod1x = [ zeros(size(q,1),1)    p.l1*sin(q(:,1))];
    rod1y = [ zeros(size(q,1),1)   -p.l1*cos(q(:,1))];
    rod2x = [ rod1x(:,2)            p.l1*sin(q(:,1)) + p.l2*sin(q(:,2))];
    rod2y = [ rod1y(:,2)           -p.l1*cos(q(:,1)) - p.l2*cos(q(:,2))];
    rod3x = [ rod1x(:,2)            p.l1*sin(q(:,1)) + p.l3*sin(q(:,3))];
    rod3y = [ rod1y(:,2)           -p.l1*cos(q(:,1)) - p.l3*cos(q(:,3))];
    pause(2);

    %animation
    figure;
    axissize = p.l1+p.l2+p.l3+5;
    for j = 1:size(t,1)    
        pause(0.01);
        plot(rod1x(j,:),rod1y(j,:),'r',rod2x(j,:),rod2y(j,:),'b',rod3x(j,:),rod3y(j,:),'g','linewidth',4)
        axis equal;
        axis([-axissize axissize -axissize axissize]);
        if j == 1 pause(2); end
    end

end