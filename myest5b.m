function[beta]=myest5b(y,X,b,u,ug,T,gamma)

% Variable setting 

bd02 = 1.06*std(y)*T^(-0.2);
[~,c] = size(X);
beta = zeros(2*c-1,length(ug));
S = zeros(2*c-1,2*c-1);
R = zeros(2*c-1,1);
    
for i = 1:length(ug)
   
    z = [X(:,1), X(:,2), X(:,3), X(:,2).*(u-ug(i)), X(:,3).*(u-ug(i))];
    eps = y - z*b;
    
    for t = 1:T
    
        we = abs(gamma - (eps(t) <= 0));
        kernel01 = 1/sqrt(2*pi)*exp(-0.5*((u(t)-ug(i))/bd02)^2);
        s = we*kernel01*z(t,:)'*z(t,:);
        S = s + S;
        r = we*kernel01*z(t,:)'*y(t);
        R = r + R;
  
    end
       
    beta(:,i) = S \ R;
    
end
 