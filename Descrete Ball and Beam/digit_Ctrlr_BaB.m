clc;
clear
close all
%% System Equations
global m g r I J muu
m = 0.2;r = 0.05;
I = 0.0002;J = 2;
g = 9.81;muu = 0;
[A,B,C,D]=State_Space();
%%
h = 0.2;
G = expm(A.*h);
%syms tav
%H = int(G,tav,[0 h])*B; H = vpa(H,4);
H=(h*eye(4)+1/2*A*h^2+1/6*A^2*h^3)*B;
rank(ctrb(G,H))

% des_poles = [0.1 0.1 -0.1 -0.1];
n=2;
des_poles = 1/2^(n)*[2 2 -2 -2];
K=acker(G,H,des_poles);

T = 5;
dt = 0.0001;t=0;
k = 1;
N = floor(h/dt);
X(:,k) = [0;0;2*3.14/180;-1*3.14/180];
%% LINEAR SYSTEM
% while t<T
%     if mod(k,N)==1
%        u(k)=-K*X(:,k);
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
       u(k)=-K*X(:,k);
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
Time=0:dt:T;
% figure;
subplot(2,2,1);plot(Time,X(1,:));title('X');
xlabel('Time');ylabel('X(t)');
subplot(2,2,2);plot(Time,X(2,:));title('X-dot');
xlabel('Time');ylabel('X-dot(t)');
subplot(2,2,3);plot(Time,X(3,:));title('Teta');
xlabel('Time');ylabel('Teta(t)');
subplot(2,2,4);plot(Time,X(4,:));title('Teta-dot');
xlabel('Time');ylabel('Teta-dot(t))');
figure;plot(Time(1:end-1),u);title('Control Effort');
xlabel('Time');ylabel('u(t)');

figure;plot(Time,X(1,:));title('X');
xlabel('Time');ylabel('X(t)');
figure;plot(Time,X(2,:));title('X-dot');
xlabel('Time');ylabel('X-dot(t)');
figure;plot(Time,X(3,:));title('Teta');
xlabel('Time');ylabel('Teta(t)');
figure;plot(Time,X(4,:));title('Teta-dot');
xlabel('Time');ylabel('Teta-dot(t))');
