function[b]=myest(y,X,u,ug,T)

% Variable setting 
bd02   = 1.06*std(y)*T^(-0.2);
[~, c] = size(X);
b      = zeros(2*c-1, length(ug));
S      = zeros(2*c-1, 2*c-1);
R      = zeros(2*c-1, 1);
 
% obatain the initial value of estimators 
for i = 1:length(ug)
    z = [X(:,1), X(:,2), X(:,3), X(:,2).*(u-ug(i)), X(:,3).*(u-ug(i))];
    for t = 1:T
        kernel01 = 1/sqrt(2*pi)*exp(-0.5*((u(t)-ug(i))/bd02)^2);
        s        = kernel01*z(t,:)'*z(t,:);
        S        = s + S;
        r        = kernel01*z(t,:)'*y(t);
        R        = r + R;
    end
    b(:, i) = S\R;
end

return(b)
