function Animate_BaB(X,t)
    L = 8;
    R = L/8;
    for i=1:length(t)
        x = X(i,1);
        teta = X(i,3);
        %% BEAM
        x1 = -L/2*cos(teta);y1=-L/2*sin(teta);
        x2 = L/2*cos(teta);y2=L/2*sin(teta);
        %% BALL
        x3 = x*cos(teta)-R*sin(teta);y3 = x*sin(teta)+R*cos(teta);
        t1 = 0:0.01:2*pi;
        x33 = x3+R*cos(t1);y33 = y3+R*sin(t1);
        %% plot
        plot(0,0,'ro',[x1 x2],[y1 y2],'b-','linewidth',5)
        hold on
        x4=[-L/2 L/2];y4=zeros(1,length(x4));
        plot(x4,y4,'-.')
        plot(x3,y3,'go',x33,y33,'g-','linewidth',3)
        axis(1/2*[-L-0.5 L+0.5 -L-0.5 L+0.5])
        axis equal
        text(-1/2*L-1, 1/2*L,['Time= ' num2str(t(i)) ' (s)'])
        text(-1/2*L-1, 1/2*L-0.5,['X= ' num2str(X(i,1)) ' (m)'])
        text(-1/2*L-1, 1/2*L-1,['Theta= ' num2str(X(i,3)/3.14*180) ' (deg)'])
        pause (0.01)
        hold off
    end
end
