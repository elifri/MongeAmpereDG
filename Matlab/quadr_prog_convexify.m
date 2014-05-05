
%load condition matrix C and previous solution
%init_C;

%cost functional
J = @(c) (solution-c)'*(solution-c);

n = size(solution,1);

H = eye(n,n);
f = -solution;

x0 =x;
%solve quadratic program
[x,fval] = quadprog(H,f,full(C),zeros(size(C,1), 1));