%{
The algorithm adds the positive constraints of varibales into the objective, and starts from a
positive primal and dual solution in which KKT condition is satisfied with
certain parameter beta=0.25. Then the algorithm iteratively calculate the descent
direction by linearlize the objective function using the first three terms of its Taylor
expansion(Newton step). By graduallly decrease the value of barrier
coefficent \mu, the objective value converges to the original optimal
objective value(as well as the solution).
%}
clear
%input data
% A = [1,1,1,0; -1,1,0,1];
% b = [2;1];
% c = [-1; -2; 0; 0];

A = [1,2,1,0,0;4,0,0,1,0;0,4,0,0,1];
b = [8; 16; 12];
c = [-2; -3; 0; 0; 0];

shape = size(A);
m = shape(1);
n = shape(2);
M = 10;%artificial variable
%multiplier = (n+2) / ((m * max(max(A)) )^m*n);%standard coefficient, each variable <= (mU)^m) 
temp = b./A;
multiplier = (n+2)/(sum(min(temp)));%standard coefficient, each variable <= min(b/A_j)
%reformulation, enlarge feasible domain by 'multiplier' times
b_bar = b*multiplier;
c = [c; M; 0];
A = [A,b_bar-A*ones(n,1), zeros(m,1);ones(1,n), 1, 1];
b = [b_bar; n+2];

%initial solution, which satisfies x>0, s>0,tolerance satisfy norm(X*S*e/u-e)=beta
u = 4 * sqrt(norm(c(1:n))^2 + M^2);
x = ones(n+2, 1);
p = [zeros(m,1); -u];
s = [c(1:n)+u*ones(n,1); u+M; u];
%paramter setting
n = n+2;%two additional variable
zeta = 0.01;%tolerance
beta = 0.25;%initial gap
k = 1;%iteration counter
K = ceil((sqrt(beta) + sqrt(n)) / (sqrt(beta) - beta)...
    *log([u-1; u-2; u; u; u+M; u]'*ones(6, 1)*(1+beta) / (zeta*(1-beta))));%maximum iteration times
alpha = 1 - (sqrt(beta) - beta)/ (sqrt(beta) + sqrt(n));%step
e = ones(n,1);
while(1)
    if (s'*x < zeta)
        break
    end
    X = diag(x);
    u = alpha*u;
    x1(k) = x(1,1)/multiplier;
    x2(k) = x(2,1)/multiplier;
    d = (eye(n) - X*X*A'*(A*X*X*A')^(-1)*A) * (X*e-X*X*c/u);
    p1 = (A*X*X*A')^(-1)*A*(X*X*c-u*X*e);
    x = x + d;
    p = p1;
    s = c - A'*p1;
    k = k+1;
end
x = x/multiplier;
x_c = 0:0.01:4;
% y1 = -x_c;
% y2 = 1+x_c;
y1 = 4-0.5*x_c;
y2 = 3;
plot(x_c, y1, 'r', x_c, y2, 'g', x1, x2, 'b*')
set(gca, 'xlim', [0,5]);
set(gca, 'ylim', [0,5]);
