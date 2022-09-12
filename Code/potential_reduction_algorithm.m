%{
The potential reduction algorithm solves the LP by introducing a potential function that combine the two
conflicting objectives, which are:
    1.Iteratively decrease the objective function value
    2.Staying away from the boundary of the feasible set
The potential function is approximated by its first order Taylor expansion
and the form is like which in the affine scaling algorithm, which has a
closed form solution. Besides, the algorithm optimizes the problem from
both the primal and dual perspective in order to take a bigger reduction in
objective value.
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
n = n+2;%one additional variable
e = ones(n,1);
zeta = 0.01;%optimality tolerance
beta = 0.285;%stepsize
gamma = 0.479;
q = n+sqrt(n);%q is large than n
k = 1;%iteration counter

while (1)
    x1(k) = x(1,1);%record coordination of x-1 axis
    x2(k) = x(2,1);%record coordination of x-2 axis
    %optimality check
    if (s'*x < zeta)
        break
    end
    X = diag(x);
    A_bar = (A*X)'/(A*X*X*A')*A*X;
    u = (eye(n)-A_bar)*(X*s*q/(s'*x)-e);
    d = -beta*X*u/(u'*u);
    %primal step
    if(u'*u >= gamma)
        x = x+d;
    %dual step
    else
        s = s'*x/q*X^(-1)*(u+e);
        p = p+(A*X*X*A')\A*X*(X*s-s'*x*e/q);
    end
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


