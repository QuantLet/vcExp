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
ng = 5;                % number of initial grid point
nreg = 3;             % number of regressors(including 1!!!)
gamma = 0.25;   % probability level
TT = [200 400];          % time intervals 
RASE = zeros(M,nreg); 
RASE0 = zeros(2*length(TT),nreg);

%% monte carlo simulations for M times

parpool('local',4);             % parallel computing tools
 
for n = 1:length(TT)
     
   T = TT(n);
   
   parfor j = 1:M
%%  Data generating and Initial setting
 
        
        [y, X] = mydgp5a(T);           % data generating process
        u = X(:,4);
        X = X(:,1:3);
        U = linspace(1,2,200)';
        ui = [U(20); U(60); U(100);U(140); U(180)];  % initial grid point 
        b = myest5d(y,X,u,ui,T);                                   % get initial value of beta from gamma = 0.5 for the initial grid point
   
%%  Estimate the functional coeffient of all the grid point
       % some lines for plotting
       % bp = zeros(1,length(u0));
       % bp0 = bp;  
       
       beta = myest5c(y,X,b,u,ui,U,T,nreg,gamma);                   % use IWLLS + one-step algorithm to calculate the parameters

%%  Calculate the RASE(root average squared error)
               
           beta0 = [3*exp(-4*(U-1).^2)+2*exp(-5*(U-2).^2), sin(sqrt(2)*pi*U), cos(sqrt(2)*pi*U)]';
         
  %%     
       RASE(j,:) = (1/T*sum((beta(1:nreg,:) - beta0).^2,2)).^0.5;            % Calculate the root average squared error
    
   end
   
%% Store the parameters of each simulation
  
  RASE0(2*n-1,:) = median(RASE);
  RASE0(2*n,:) = std(RASE); 

end

% plot(u0,bp0,'k-',u0,bp,'k.')

delete(gcp)