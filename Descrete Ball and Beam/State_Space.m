function [A,B,C,D]=State_Space()
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

    A = vpa(A,6);B = vpa(B,6);
    A = double(A);B = double(B);
    C = [1 0 0 0];D = 0;
end

