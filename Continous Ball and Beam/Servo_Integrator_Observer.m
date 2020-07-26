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
%% DESIGN STATE_FEEDBACK CONTROLLER
dt = 0.01;
T = 50;
Abar = [A zeros(4,1);-C 0];
Bbar = [B;0];
BB= [zeros(4,1);1];
CC = [C 0];
%% DESIGN STATE FEEDBACK CONTROLLER
desired_poles = [-2+1j -2-1j -2 -2 -4];
K = acker(Abar,Bbar,desired_poles);
Acl = Abar-Bbar*K;
%% DESIGN STATE_OBSERVER
rank(obsv(A,C))
Desired_poles = [-5 -5 -5 -5];
L = acker(A',C',Desired_poles)';

%% SOLVE EQUATIONS
in = input("Inter 0 or 1:");
init = [0 0 -5*3.14/180 2*3.14/180 0 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:dt:T;
if in==0
    [t,X] = ode45(@(t,x) linear_ode(t,x,A,B,Abar,Bbar,BB,K,L),tspan,init,options);
else 
    [t,X] = ode45(@(t,x) nonlinear_ode(t,x,K,L),tspan,init,options);  
end
%% PLOTs
plot(t,X(:,1),'b',t,X(:,5)/2.04,'r--',t,0.3*sign(sin(0.5*t)),'g')
title('X')
xlabel('Time');ylabel('X')
%% FUNCTIONS
function dx = nonlinear_ode(t,XX,K,L)
    global m r g I J mu
    Yr = 0.3*sign(sin(0.5*t));
    X = XX(1:5);
    Xhat=XX(6:9);
    u = -K*X;
    
    f1=1/(J+I+m*(X(1)^2+r^2))*(m*g*(cos(X(3))*X(1)+sin(X(3))*r)-u);
    f2=m*g*cos(X(3))-m*(X(1)*f1+2*X(2)*X(4));
    dX=[X(2);
        X(1)*X(4)^2+g*sin(X(3))-mu/m*f2;
        X(4);
        f1;
        Yr-X(1)];
    
    f1=1/(J+I+m*(Xhat(1)^2+r^2))*(m*g*(cos(Xhat(3))*Xhat(1)+sin(Xhat(3))*r)-u);
    f2=m*g*cos(Xhat(3))-m*(Xhat(1)*f1+2*Xhat(2)*Xhat(4));
    dXhat=[ Xhat(2);
            Xhat(1)*Xhat(4)^2+g*sin(Xhat(3))-mu/m*f2;
            Xhat(4);
            f1];
        
    dXhat = dXhat + L*(Yr-Xhat(1));
    dx = [dX;dXhat];
end
function dx = linear_ode(t,XX,A,B,Abar,Bbar,BB,K,L)
    Yr = 0.3*sign(sin(0.5*t));
    
    X = XX(1:5);
    Xhat=XX(6:9);
    u = -K*X;
    dX = Abar*X + Bbar*u + [zeros(4,1);Yr];
    Yhat = Xhat(1);
    dXhat = A*Xhat + B*u + L*(Yr-Yhat);
    dx=[dX;dXhat];
end


