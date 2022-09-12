%{
The affine scaling algorithm start from a primal positive solution, and
iterately calculate the descent direction by solving a SOCP with the same
objective function as original formulations within a ellipsoid, where the SOCP
has a closed form solution.
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
M = 10;%a big number
%reformulation and initialization to obtain a positive solution
A = [A,b-A*ones(n,1)];
c = [c; M];
x = [ones(n, 1); 1];

%parameters of the algorithm
n = n+1;%one additional variable
e = ones(n,1);
zeta = 0.01;%optimality tolerance
beta = 0.5;%stepsize
k = 1;%iteration counter

while(1)
    x1(k) = x(1,1);%record coordination of x-1 axis
    x2(k) = x(2,1);%record coordination of x-2 axis
    %computation of dual estimates and reduced costs
    X = diag(x);
    p = (A*X*X*A')\A*X*X*c;%dual estimates
    r = c-A'*p;%reduced costs
    %optimality check
    if (all(r >= 0) && e'*X*r<zeta)
        break;
    end
    %unboundedness check
    if (all(-X*X*r >= 0))
        disp("Problem unbounded")
        break
    end
    %x = x - beta*X*X*r/norm(X*r);%standard stepsize
    x = x - beta*X*X*r/max(abs(X*r));
    %x = x - beta*X*X*r/max(X*r);
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



