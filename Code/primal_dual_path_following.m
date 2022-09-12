%{
The algorithm starts from a feasible solution which both of the primal variables and dual
slack variables are positively by modifying the original formulation with two auxiliary 
variables introduced.
The code solves a linear programming problem(LP) using primal-dual path following algorithm, 
which tries to solve the system of equations defined by KKT optimality
condition using Newton's method.
The code accepts an input of coefficient matrix A, resouce constraints vector b, 
and cost vector c.
%}

clear
%input data
% A = [1,1,1,0; -1,1,0,1];
% b = [2;1];
% c = [-1; -2; 0; 0];

A = [1,2,1,0,0;4,0,0,1,0;0,4,0,0,1];
b = [8;16;12];
c = [-2;-3;0;0;0];

shape = size(A);
m = shape(1);
n = shape(2);
%reformulation and initial solutions, 2 variables of primal problem are introdcued
%initial solution, with x>0, s>0
M1 = 30;%artificial variable
M2 = 30;
e = ones(n,1);
x = [ones(n, 1); 1; M2-(e-c)'*e];
p = [zeros(m,1); -1];
s = [e; M1; 1];
A = [A,b-A*ones(n,1), zeros(m,1);ones(1,n)-c', 0, 1];
b = [b; M2];
c = [c; M1; 0];

%parameters of the algorithm
n = n+2;%two additional variable
zeta = 0.01;%tolerance
rho = 0.7;%decide barrier parameter by 'rho*(s'*x)/n'
e = ones(n,1);
alpha = 0.9;%avoid the point moving to the boundary
k = 1;%iteration counter
while(1)
    if (s'*x < zeta)
        break
    end
    %disp(s'*x);
    x1(k) = x(1,1);%record coordination of x-1 axis
    x2(k) = x(2,1);%record coordination of x-2 axis
    u = rho * (s'*x) / n;%%barrier parameter
    X = diag(x);
    S = diag(s);
    %solve equation system
    D = sqrt(X/S);
    P = D*A'/(A*D*D*A')*A*D;
    vu = X\D*(u*e-X*S*e);
    %descent direction
    d_x = D*(eye(n)-P)*vu;
    d_p = -(A*D*D*A')\A*D*vu;
    d_s = D\P*vu;
    %find step lengths
    negativeP = x./d_x;
    betaP = min(1, -max(negativeP(negativeP<0))*alpha);
    negativeD = x./d_s;
    betaD = min(1, -max(negativeD(negativeD<0))*alpha);
    %update solution vectors
    x = x + betaP*d_x;
    p = p + betaD*d_p;
    s = s + betaD*d_s;
    k = k+1;
end
%show convergence curve
x_c = 0:0.01:4;
% y1 = -x_c;
% y2 = 1+x_c;
y1 = 4-0.5*x_c;
y2 = 3;
plot(x_c, y1, 'r', x_c, y2, 'g', x1, x2, 'b*')
set(gca, 'xlim', [0,5]);
set(gca, 'ylim', [0,5]);