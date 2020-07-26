clc;
clear
close all
%% System Equations
global m g r I J muu
m = 0.2;r = 0.05;
I = 0.0002;J = 2;
g = 9.81;muu = 0;
[A,B,C,D]=ball_and_beam_ss();
h = 0.1;
G = expm(A*h);
%syms tav
% H = int(G,tav,[0 h])*B; H = vpa(H,4);
H=(h*eye(4)+1/2*A*h^2+1/6*A^2*h^3+1/24*A^3*h^4+1/120*A^4*h^5)*B;
rank(ctrb(G,H))

des_poles = [0.5 0.5 -0.5 -0.5];
K=acker(G,H,des_poles);

rank(obsv(G,C))
des_poles = [0.1 0.1 -0.1 -0.1];
L=acker(G',C',des_poles)';

T = 10;
dt = 0.0001;
t=0;
k = 1;k1 = 0;
N = floor(h/dt);
X(:,k) = [0;0;2*3.14/180;-1*3.14/180];
Xh(:,k) = [0;0;0;0];
%% LINEAR SYSTEM
% while t<T
%     if mod(k,N)==1
%        k1 = k1+1;
%        u(k)=-K*Xh(:,k1);
%        Y=C*X(:,k);
%        Yh=C*Xh(:,k1);
%        Xh(:,k1+1)=G*Xh(:,k1)+H*u(k)+L*(Y-Yh);  
%     else
%        u(k)=u(k-1);
%     end
%     X(:,k+1)=X(:,k)+(A*X(:,k)+B*u(k))*dt;
%     k = k+1;
%     t=t+dt;    
% end
%% NONLINEAR SYSTEM
while t<T
    if mod(k,N)==1
       k1 = k1+1;
       u(k)=-K*Xh(:,k1);
       Y=C*X(:,k);
       Yh=C*Xh(:,k1);
       Xh(:,k1+1)=G*Xh(:,k1)+H*u(k)+L*(Y-Yh);  
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
    k = k+1;
    t=t+dt;    
end
%% PLOTS
Time=0:dt:T+dt;
Time1=0:h:T+h;
plot(Time,X(1,:),Time1,Xh(3,:));title('X');
xlabel('Time');ylabel('X(t)');
legend('WITHOUT observer','WITH observer');
figure;plot(Time,X(3,:),Time1,Xh(3,:));title('Teta');
xlabel('Time');ylabel('Teta(t)');
legend('WITHOUT observer','WITH observer');
figure;plot(Time(1:end-1),u);title('Control Effort');
xlabel('Time');ylabel('u(t)');
%% FUNCTIONS
function [A,B,C,D]=ball_and_beam_ss()
    syms x1 x2 x3 x4 u
    global m g r I J muu
    f1=1/(J+I+m*(x1^2+r^2))*(m*g*(cos(x3)*x1+sin(x3)*r)-u);
    f2=m*g*cos(x3)-m*(x1*f1+2*x2*x4);
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