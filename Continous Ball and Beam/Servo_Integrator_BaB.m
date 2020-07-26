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
mu = 0;
[A,B,C,D]=State_Space();
% DESIGN STATE_FEEDBACK CONTROLLER
dt = 0.01;
T = 50;
AA = [A zeros(4,1);-C 0];
BB = [B;0];
BBB= [zeros(4,1);1];
CC = [C 0];
% DESIGN STATE FEEDBACK CONTROLLER
desired_poles = [-2+1j -2-1j -2 -2 -4];
KK = acker(AA,BB,desired_poles);
Acl = AA-BB*KK;
%% Plots
in = input("Inter 0 or 1:");
init = [0.1 0.5 -5*3.14/180 2*3.14/180 0];
% init = [0.1 0.5 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:dt:T;
if in==0
    [t,X] = ode45(@(t,x) linear_ode(t,x,AA,BB,BBB,KK),tspan,init,options);
else 
    [t,X] = ode45(@(t,x) nonlinear_ode(t,x,KK),tspan,init,options);  
end
plot(t,X(:,1),t,0.3*sign(sin(0.5*t)))
title('X')
xlabel('Time');ylabel('X')
%%
function dX = nonlinear_ode(t,X,K)
    global m r g I J mu
    Yr = 0.3*sign(sin(0.5*t));
    u = -K*X;
    
    f1=1/(J+I+m*(X(1)^2+r^2))*(m*g*(cos(X(3))*X(1)+sin(X(3))*r)-u);
    f2=m*g*cos(X(3))-m*(X(1)*f1+2*X(2)*X(4));
    dX=[X(2);
        X(1)*X(4)^2+g*sin(X(3))-mu/m*f2;
        X(4);
        f1;
        Yr-X(1)];
end
function dX = linear_ode(t,X,A,B,BB,K)
    Yr = 0.3*sign(sin(0.5*t));
    u = -K*X;
    dX = A*X + B*u + BB*Yr;
end
