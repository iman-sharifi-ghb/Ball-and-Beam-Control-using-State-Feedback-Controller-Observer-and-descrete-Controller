clc;
clear
close all
%% System Equations
global m g r I J mu
m = 0.2;
r = 0.05;
g = 9.81;
I = 0.0002;
J = 2;
mu=0;
[A,B,C,D]=State_Space();
% DESIGN STATE_FEEDBACK CONTROLLER
rank(ctrb(A,B))
desired_poles = [-2+1j -2-1j -2 -2];
K = acker(A,B,desired_poles);
%% LINEAR ODE45
% init = [0.1 0.5 -5 2];
init = [1 0.5 pi/10 -0.5];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:0.01:10;
[t,X] = ode45(@(t,x) linear_ode(t,x,A,B,K),tspan,init,options);
% subplot(2,2,1);plot(t,X(:,1));title('X');
% subplot(2,2,2);plot(t,X(:,2));title('X-dot');
% subplot(2,2,3);plot(t,X(:,3)/3.14*180);title('Teta');
% subplot(2,2,4);plot(t,X(:,4)/3.14*180);title('Teta-dot');
% figure;
% [t,XX] = ode45(@(t,x) nonlinear_ode(t,x,K),tspan,init,options);
% subplot(2,2,1);plot(t,XX(:,1));title('X');
% subplot(2,2,2);plot(t,XX(:,2));title('X-dot');
% subplot(2,2,3);plot(t,XX(:,3)/3.14*180);title('Teta');
% subplot(2,2,4);plot(t,XX(:,4)/3.14*180);title('Teta-dot');
% figure;
% plot(t,X(:,1),t,XX(:,1),'r')
% legend('Linear Equation','NonLinear Equation')
% title('Compare X')
% xlabel('Time');ylabel('X')
% figure;
% plot(t,X(:,3),t,XX(:,3),'r')
% legend('Linear Equation','NonLinear Equation')
% title('Compare Teta')
% xlabel('Time');ylabel('Teta')
%% ANIMATION
t=0:0.01:5;
Animate_BaB(X,t)
%%
function dx = nonlinear_ode(t,x,K)
    global m r g I J mu
    u = -K*x;
    f1=1/(J+I+m*(x(1)^2+r^2))*(m*g*(cos(x(3))*x(1)+sin(x(3))*r)-u);
    f2=m*g*cos(x(3))-m*(x(1)*f1+2*x(2)*x(4));
    dx=[x(2);
        x(1)*x(4)^2+g*sin(x(3))-mu/m*f2;
        x(4);
        f1];
end
function dx = linear_ode(t,x,A,B,K)
    u = -K*x;
    dx = A*x + B*u;
end