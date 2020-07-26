clc;
clear;close all
%   X(k+1)=GX(k)+Hu(k)
%    y(k)=CX(k)+Du(k)
%%  Servo Type 1 + Observer
global m g r I J muu
m = 0.2;r = 0.05;
I = 0.0002;J = 2;
g = 9.81;muu = 0;

[A,B,C,D]=ball_and_beam_ss();
h = 0.2;
G = expm(A*h);
% syms tav
% H = int(G,tav,[0 h])*B; H = vpa(H,4);
H=(h*eye(4)+1/2*A*h^2+1/6*A^2*h^3+1/24*A^3*h^4+1/120*A^4*h^5)*B;
rank(ctrb(G,H))
rank(obsv(G,C))

rank(obsv(G,C))
des_poles = [0.1 0.1 -0.1 -0.1];
L=acker(G',C',des_poles)';

mu=[0.1 0.1 -0.1 -0.1];
K=acker(G,H,mu);
eig(G-H*K)

Gcl0=-C/(G-H*K)*H;
X(:,1)=[0.5;0;5*pi/180;-1*pi/180];
N=100;
dt = 0.0001;
%% DIGITAL SYSTEM
for k=1:N-1
    yd(k)=2*sign(sin(0.2*k));
    u(k)=-K*X(:,k)+yd(k)/Gcl0;
    X(:,k+1)=G*X(:,k)+H*u(k);
    y(k)=C*X(:,k);
    err=yd(k)-y(k);
    psi(k+1)=psi(k)+err;
end
%% LINEAR SYSTEM
% for k=1:N-1
%     yd(k)=2*sign(sin(0.2*k));
%     u(k)=-K*X(:,k)+yd(k)/Gcl0;
%     X(:,k+1)=X(:,k)+(A*X(:,k)+B*u(k))*dt;
%     y(k)=C*X(:,k);
%     err=yd(k)-y(k);
%     psi(k+1)=psi(k)+err;
% end
%% NONLINEAR SYSTEM
for k=1:N-1
    yd(k)=2*sign(sin(0.2*k));
    u(k)=-K*X(:,k)+yd(k)/Gcl0;
    f1=1/(J+I+m*(X(1,k)^2+r^2))*(m*g*(cos(X(3,k))*X(1,k)+sin(X(3,k))*r)-u(k));
    f2=m*g*cos(X(3,k))-m*(X(1,k)*f1+2*X(2,k)*X(4,k));
    SYS = [X(2,k);
           X(1,k)*X(4,k)^2+g*sin(X(3,k))-muu/m*f2;
           X(4,k);
           f1];
    X(:,k+1)=X(:,k)+(SYS)*dt;
    y(k)=C*X(:,k);
    err=yd(k)-y(k);
    psi(k+1)=psi(k)+err;
end
%% PLOTS
Time=1:N;
X(1,:)=4.*X(1,:).*sin(0.2*Time);
plot(Time,X(1,:),'b',Time(1:N-1),yd,'r');title('X');
xlabel('Time');ylabel('X(t)');
figure;
plot(u);title('Control Effort');
xlabel('Time');ylabel('u(t)');
%% ANIMATION
% t=0:0.01:5;
% animate_ball_beam(X,t)
%% FUNCTIONS
function [A,B,C,D]=ball_and_beam_ss()
    syms x1 x2 x3 x4 u
    global m g r I J muu
    f1=1/(J+I+m*(x1^2+r^2))*(m*g*(cos(x3)*x1+sin(x3)*r)-u);
    f2 = m*g*cos(x3)-m*(x1*f1+2*x2*x4);
    f1 = simplify(f1);f2 = simplify(f2);
    dx1 = x2;
    dx2 = x1*x4^2+g*sin(x3)-muu/m*f2;
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

    A = vpa(A,6);B = vpa(B,6);
    A = double(A);B = double(B);
    C = [1 0 0 0];D = 0;
end