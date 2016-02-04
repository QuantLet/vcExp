function[beta]=myest(y, X, b, u, ug, T, gamma)

% Variable setting 
bd01   = 1.06*std(y)*T^(-0.2);
[~, c] = size(X);
S      = zeros(2*c - 1, 2*c - 1);
R      = zeros(2*c - 1, 1);

% Calculate the initial estimators
bx = b;
b0 = zeros(2*c - 1,1);
z  = [X(:, 1),X(:, 2),X(:, 3), X(:, 2).*(u - ug), X(:, 3).*(u - ug)];
k  = 0;

% Iterate to obatain the estimator of convergence
while((norm(abs(b0 - bx)) > 1e-5)&&(k <= 1000))   
    eps = y - z*bx;
    for t = 1:T
        we       = abs(gamma - (eps(t) <= 0));
        kernel01 = 1/sqrt(2*pi)*exp(-0.5*((u(t) - ug)/bd01)^2);
        s        = we*kernel01*z(t, :)'*z(t, :);
        S        = s + S;
        r        = we*kernel01*z(t, :)'*y(t);
        R        = r + R;
    end
   bs = S\R;
   b0 = bx;
   bx = bs;
   k  = k+1;
end

% Calculate the final resultx
beta = bx;
return(beta)
