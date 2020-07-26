clc;
clear
close all
%% System Equations
global M m l g mu
M = 5;
m = 1;
l = 0.5;
g = 9.81;
mu = 0;
[A,B,C,D]=State_Space();
%%
dt = 0.01;
T = 50;
t = 0:dt:T;
Yr = 0.3*sign(sin(0.5*t));
% DESIGN STATE FEEDBACK CONTROLLER
desired_poles = [-2+1j -2-1j -5 -5];
K = acker(A,B,desired_poles);
Acl = A-B*K;
G0 = -C/(A-B*K)*B;
%% Plots
in = input("Inter 0 or 1:");
init = [0.1 0.5 -5*3.14/180 2*3.14/180];
% init = [0.1 0.5 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:dt:T;
if in==0
    [t,X] = ode45(@(t,x) linear_ode(t,x,A,B,C,K),tspan,init,options);
else 
    [t,X] = ode45(@(t,x) nonlinear_ode(t,x,G0,K),tspan,init,options);  
end
plot(t,X(:,1),t,0.3*sign(sin(0.5*t)))
title('X')
xlabel('Time');ylabel('X')

function dX = linear_ode(t,X,A,B,C,K)
    Yr = 0.3*sign(sin(0.5*t));
    G0 = -C/(A-B*K)*B;
    u = -K*X+Yr/G0;
    dX = A*X + B*u;
end
function dX = nonlinear_ode(t,X,G0,K)
    global m r g I J mu
    Yr = 0.3*sign(sin(0.5*t));
    u = -K*X+Yr/G0;
    
    f1=1/(J+I+m*(X(1)^2+r^2))*(m*g*(cos(X(3))*X(1)+sin(X(3))*r)-u);
    f2=m*g*cos(X(3))-m*(X(1)*f1+2*X(2)*X(4));
    dX=[X(2);
        X(1)*X(4)^2+g*sin(X(3))-mu/m*f2;
        X(4);
        f1];
end
