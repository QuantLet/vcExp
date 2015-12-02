% Title: A Semiparametric Estimation of Patially Varying-Coefficient
% Expectile Model
% Date: Dec 15, 2015
% Author: Zongwu Cai, Ying Fang and Dingshi Tian

% clearing work&preparing
clc
clear
close all

% parameter setting

M = 100;             % number of simulations
C = 1;                  % choose the different DGP
ng = 5;                % number of initial grid point
nreg = 3;             % number of regressors(including 1!!!)
gamma = 0.25;   % probability level
TT = [200];          % time intervals 
RASE = zeros(M,nreg); 
RASE0 = zeros(2*length(TT),nreg);

%% monte carlo simulations for M times

parpool('local',4);             % parallel computing tools
 
for n = 1:length(TT)
     
   T = TT(n);
   
   parfor j = 1:M
%%  Data generating and Initial setting
 
        
        [y, X] = mydgp5(T,C);           % data generating process
       
        if C == 1   
            u = X(:,2);                          % smooth variable
        elseif C ==2
            u = X(:,3);
        else 
            u = X(:,4);
            X = X(:,1:3);
        end
        
        b = myest5d(y,X,u,u,T);  
        beta = zeros(2*nreg-1, length(u));                     
   
%%  Estimate the functional coeffient of all the grid point
       % some lines for plotting
       % bp = zeros(1,length(u0));
       % bp0 = bp;  
      
       for i = 1:length(u)
        
            beta(:,i) = myest5a(y,X,b(:,i),u,u(i),T,gamma);  % use IWLLS to calculate the parameters for initial point
        
       end
                      % use IWLLS + one-step algorithm to calculate the parameters

%%  Calculate the RASE(root average squared error)
   
       if C == 1
            
           beta0 = [zeros(T,1), (1-0.18*(u>0)).*(0.138 + (0.316+0.982*u).*exp(-3.89*u.^2)), -0.437 - (0.659+0.126*u).*exp(-3.89*u.^2)]';  % real beta   
       
       elseif C == 2
           
           beta0 = [zeros(T,1), 0.4*(U <= 1) - 0.8*(U > 1), -0.6*(U <= 1) + 0.8*(U > 1)]';
           
       else
           
           beta0 = [zeros(T,1), sin(sqrt(2)*pi*U), cos(sqrt(2)*pi*U)]';
       
       end
  %%     
       RASE(j,:) = (1/T*sum((beta(1:nreg,:) - beta0).^2,2)).^0.5;            % Calculate the root average squared error
    
   end
   
%% Store the parameters of each simulation
  
  RASE0(2*n-1,:) = median(RASE);
  RASE0(2*n,:) = std(RASE); 

end

% plot(u0,bp0,'k-',u0,bp,'k.')

delete(gcp)
