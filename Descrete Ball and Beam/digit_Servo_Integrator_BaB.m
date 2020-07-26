clc;
clear;close all
%   X(k+1)=GX(k)+Hu(k)
%    y(k)=CX(k)+Du(k)
%%  Servo Type 1 + Observer
global m g r I J muu
m = 0.2;r = 0.05;
I = 0.0002;J = 2;
g = 9.81;muu = 0;

[A,B,C,D]=State_Space();
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

mu=0*[0.1 0.1 -0.1 -0.1 0];
Gb=[G zeros(4,1);-C 1];
Hb=[H;0];
rank(ctrb(Gb,Hb))

Kb=acker(Gb,Hb,mu);
eig(Gb-Hb*Kb)
X(:,1)=0.1*[0.5;0;5*pi/180;-1*pi/180];
psi(1)=0;
N=100;
t=0;
dt = 0.0001;
%% DIGITAL SYSTEM
% for k=1:N-1
%     yd(k)=2*sign(sin(0.2*k));
%     u(k)=-Kb(1:4)*X(:,k)-Kb(5)*psi(k);
%     X(:,k+1)=G*X(:,k)+H*u(k);
%     y(k)=C*X(:,k);
%     err=yd(k)-y(k);
%     psi(k+1)=psi(k)+err;
% end
%% LINEAR SYSTEM
% for k=1:N-1
% if 
%     yd(k)=2*sign(sin(0.2*k));
%     u(k)=-Kb(1:4)*X(:,k)-Kb(5)*psi(k);
%     X(:,k+1)=X(:,k)+(A*X(:,k)+B*u(k))*dt;
%     y(k)=C*X(:,k);
%     err=yd(k)-y(k);
%     psi(k+1)=psi(k)+err;
% end
%% LINEAR SYSTEM
% for k=1:N-1
%     yd(k)=0.1*sign(sin(0.1*k));
%     if mod(k,N)==1
%        u(k)=-Kb(1:4)*X(:,k)-Kb(5)*psi(k);
%     else
%        u(k)=u(k-1);
%     end
%     X(:,k+1)=X(:,k)+(A*X(:,k)+B*u(k))*dt;
%     y(k)=C*X(:,k);
%     disp(A*X(:,k)+B*u(k))
%     err=yd(k)-y(k);
%     psi(k+1)=psi(k)+err;
% end
%% NONLINEAR SYSTEM
for k=1:N-1
    yd(k)=0.2*sign(sin(0.1*k));
    if mod(k,N)==1
       u(k)=-Kb(1:4)*X(:,k)-Kb(5)*psi(k);
    else
       u(k)=u(k-1);
    end
    f1=1/(J+I+m*(X(1,k)^2+r^2))*(m*g*(cos(X(3,k))*X(1,k)+sin(X(3,k))*r)-u(k));
    f2=m*g*cos(X(3,k))-m*(X(1,k)*f1+2*X(2,k)*X(4,k));
    SYS = [X(2,k);
           X(1,k)*X(4,k)^2+g*sin(X(3,k))-muu/m*f2;
           X(4,k);
           f1];
    X(:,k+1)=X(:,k)+(SYS)*dt;
    y(k)=C*X(:,k);
    err=yd(k)-y(k);
    psi(k+1)=psi(k)+err*dt;
    k = k + 1;
end
%% PLOTS
Time=1:N;
% X(1,:)=4.*X(1,:).*sin(0.2*Time);
plot(Time,X(1,:),'b',Time(1:N-1),yd,'r');title('X');
xlabel('Time');ylabel('X(t)');
figure;
plot(u);title('Control Effort');
xlabel('Time');ylabel('u(t)');
%% ANIMATION
% t=0:0.01:5;
% animate_ball_beam(X,t)
