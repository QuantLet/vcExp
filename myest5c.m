function [beta] = myest5c(y,X,b,u,ui,U,T,nreg,gamma)
%%
 ng = 5; 
 beta = zeros(2*nreg-1, length(U));
 lg = length(U);
   
 %% estimate the parameter of some initial point
        
for i = 1:ng
        
    beta(:,(i-0.5)*lg/ng) = myest5a(y,X,b(:,i),u,ui(i),T,gamma);  % use IWLLS to calculate the parameters for initial point
        
end

%% use the one step algorithm to caculate the remaining grid point       

for i = 1:ng
               
   beta(:,[(i-0.5)*lg/ng-1,(i-0.5)*lg/ng+1]) = myest5b(y,X,beta(:,(i-0.5)*lg/ng),u,...
                                                                                U([(i-0.5)*lg/ng-1,(i-0.5)*lg/ng+1]),T,gamma);            % grid point +_1 
end
        
for i = 1 : ng
            
     for k = 2 : lg/(2*ng)-1
   
          beta(:,(i-0.5)*lg/ng-k) = myest5b(y,X,beta(:,(i-0.5)*lg/ng-k+1),u,U((i-0.5)*lg/ng-k),T,gamma);          % all other grid point -
          beta(:,(i-0.5)*lg/ng+k) = myest5b(y,X,beta(:,(i-0.5)*lg/ng+k-1),u,U((i-0.5)*lg/ng+k),T,gamma);         % all other grid point +
             
     end
     
     beta(:,40*i) = myest5b(y,X,beta(:,40*i-1),u,U(40*i),T,gamma);

end