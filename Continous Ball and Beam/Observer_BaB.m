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
mu= 0;
[A,B,C,D]=disk_on_rod_ss();
% DESIGN STATE_FEEDBACK CONTROLLER
rank(ctrb(A,B))
desired_poles = [-2+1j -2-1j -2 -2];
K = acker(A,B,desired_poles);
% DESIGN STATE_OBSERVER
rank(obsv(A,C))
Desired_poles = [-5 -5 -5 -5];
L = acker(A',C',Desired_poles)';
%% DESIGN REDUCE-ORDER OBSERVER
Aaa = A(1,1);Aab = A(1,2:4);
Aba = A(2:4,1);Abb = A(2:4,2:4);
Ba = B(1);Bb = B(2:4);
rank(obsv(Abb,Aab))
Des_poles = [-5 -5 -5];
Lr = acker(Abb',Aab',Des_poles)';
%% Plots
dt = 0.01;
T = 10;
in = input("Inter 0 or 1 or 2 or 3:");
init = [0.1 0.5 -5*3.14/180 2*3.14/180 0 0 0 0];
% init = [0.1 0.5 0 0 0 0 0 0];
options = odeset('RelTol',1e-2,'AbsTol',1e-4);
tspan = 0:dt:T;
if in==0
    [t,XX] = ode45(@(t,x) linear_ode(t,x,A,B,C,K,L),tspan,init,options);
elseif in==1
    [t,XX] = ode45(@(t,x) nonlinear_ode(t,x,C,K,L),tspan,init,options);
elseif in==2
    [t,XX] = ode45(@(t,x) linear_ode_order_reduction(t,x,A,B,K,Lr),tspan,init,options);
else
    [t,XX] = ode45(@(t,x) nonlinear_ode_reduce_order(t,x,K,Lr),tspan,init,options);
end
X = XX(:,1:4);
subplot(2,2,1);plot(t,X(:,1));title('X');
subplot(2,2,2);plot(t,X(:,2));title('X-dot');
subplot(2,2,3);plot(t,X(:,3)/3.14*180);title('Teta');
subplot(2,2,4);plot(t,X(:,4)/3.14*180);title('Teta-dot');
figure
Xhat = XX(:,5:8);
subplot(2,2,1);plot(t,Xhat(:,1));title('X');
subplot(2,2,2);plot(t,Xhat(:,2));title('X-dot');
subplot(2,2,3);plot(t,Xhat(:,3)/3.14*180);title('Teta');
subplot(2,2,4);plot(t,Xhat(:,4)/3.14*180);title('Teta-dot');
figure
plot(t,X(:,1),t,Xhat(:,1))
legend('state','observer')
title('X')
xlabel('Time');ylabel('X')
figure
plot(t,X(:,2),t,Xhat(:,2))
legend('state','observer')
title('X-dot')
xlabel('Time');ylabel('X-dot')
figure
plot(t,X(:,3)/180*3.14,t,Xhat(:,3)/180*3.14)
legend('state','observer')
title('Teta')
xlabel('Time');ylabel('Teta')
figure
plot(t,X(:,4)/180*3.14,t,Xhat(:,4)/180*3.14)
legend('state','observer')
title('X4')
xlabel('Time');ylabel('X4')
%%
function dx = nonlinear_ode(t,XX,C,K,L)
    global m r g I J mu
    X = XX(1:4);
    Xhat=XX(5:8);
    u = -K*Xhat;
    Y = C*X;
    Yhat = C*Xhat;
    
    f1=1/(J+I+m*(X(1)^2+r^2))*(m*g*(cos(X(3))*X(1)+sin(X(3))*r)-u);
    f2=m*g*cos(X(3))-m*(X(1)*f1+2*X(2)*X(4));
    dX=[X(2);
        X(1)*X(4)^2+g*sin(X(3))-mu/m*f2;
        X(4);
        f1];
    
    f1=1/(J+I+m*(Xhat(1)+r)^2)*(m*g*(cos(Xhat(3))*Xhat(1)+sin(Xhat(3))*r)-u);
    f2=m*g*cos(Xhat(3))-m*(Xhat(1)*f1+2*Xhat(2)*Xhat(4));
    dXhat=[Xhat(2);
        Xhat(1)*Xhat(4)^2+g*sin(Xhat(3))-mu/m*f2;
        Xhat(4);
        f1];
    dXhat = dXhat + L*(Y-Yhat);
    dx = [dX;dXhat];
end
function dx = linear_ode(t,XX,A,B,C,K,L)
    X = XX(1:4);
    Xhat=XX(5:8);
    u = -K*Xhat;
    dX = A*X + B*u;
    Y = C*X;
    Yhat = C*Xhat;
    dXhat = A*Xhat + B*u + L*(Y-Yhat);
    dx=[dX;dXhat];
end
function [A,B,C,D]=disk_on_rod_ss()
    syms x1 x2 x3 x4 u
    global m g r I J mu
    f1=1/(J+I+m*(x1^2+r^2))*(m*g*(cos(x3)*x1+sin(x3)*r)-u);
    f2=m*g*cos(x3)-m*(x1*f1+2*x2*x4);
    dx1 = x2;
    dx2 = x1*x4^2+g*sin(x3)-mu/m*f2;
    dx3 = x4;
    dx4 = f1;
    
    x = [x1;x2;x3;x4];
    dx = [dx1;dx2;dx3;dx4];

    A = jacobian(dx,x);
    A = simplify(A);
    B = jacobian(dx,u);
    B = simplify(B);

    A = subs(A,[x1,x2,x3,x4,u],[0,0,0,0,0]);
    B = subs(B,[x1,x2,x3,x4,u],[0,0,0,0,0]);

    A = vpa(A,6);
    B = vpa(B,6);
    A = double(A);
    B = double(B);
    C = [1 0 0 0];
    D = 0;
end
function dx = linear_ode_order_reduction(t,XX,A,B,K,Lr)
    Aaa = A(1,1);Aab = A(1,2:4);
    Aba = A(2:4,1);Abb = A(2:4,2:4);
    Ba = B(1);Bb = B(2:4);
    
    Xa = XX(1);
    Xb = XX(2:4);
    Xhat_b = XX(6:8);
    u = -K*[Xa;Xhat_b];
    
    dXa = Aaa*Xa + Aab*Xb + Ba*u;
    dXb = Aba*Xa + Abb*Xb + Bb*u;
    dXhat_a = dXa;
    dXhat_b = Aba*Xa + Abb*Xhat_b + Bb*u + Lr*(Aab*Xb - Aab*Xhat_b);

    dx=[dXa;dXb;dXhat_a;dXhat_b];
end
function dx = nonlinear_ode_reduce_order(t,XX,K,Lr)
    global m r g I J mu
    Xa = XX(1);
    Xb = XX(2:4);
    Xhat_b=XX(6:8);
    u = -K*[Xa;Xhat_b];
    
    f1=1/(J+I+m*(Xa^2+r^2))*(m*g*(cos(Xb(2))*Xa+sin(Xb(2))*r)-u);
    f2=m*g*cos(Xb(2))-m*(Xa*f1+2*Xb(1)*Xb(3));
    dXa=Xb(1);
    dXb=[Xa*Xb(3)^2+g*sin(Xb(2))-mu/m*f2;
        Xb(3);
        f1];
    
    f1=1/(J+I+m*(Xa^2+r^2))*(m*g*(cos(Xhat_b(2))*Xa+sin(Xhat_b(2))*r)-u);
    f2=m*g*cos(Xhat_b(2))-m*(Xa*f1+2*Xhat_b(1)*Xhat_b(3));
    dXhat_a = dXa;
    dXhat_b=[Xa*Xhat_b(3)^2+g*sin(Xhat_b(2))-mu/m*f2;
        Xhat_b(3);
        f1];
    dXhat_b = dXhat_b + Lr*(Xb(1)-Xhat_b(1));
    dx = [dXa;dXb;dXhat_a;dXhat_b];
end







