
%load condition matrix C and previous solution
%init_C;

%cost functional
%J = @(c) (solution-c)'*(solution-c);

%n = size(solution,1);
n = size(C_values, 1);

% to minimise distance in the coefficients
%H = eye(n,n);
%f = -solution;

% to minimum distance to data points
H = A'*A;
f = -A'*C_values;

opts = optimset('Algorithm','interior-point-convex','Display','iter', 'TolX', 1e-10, 'TolCon', 1e-10);

%solve quadratic program
[x,fval] = quadprog(H,f,full(-C),zeros(size(C,1),1), [],[],[],[],[], opts);