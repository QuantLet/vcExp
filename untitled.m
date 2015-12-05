% This is an bivariate dimensional example of how to utilize Thorp Critrion to find the optimal
% fraction of an portfolio. We use Dax30 and SP500 from 12.31.2004 to 11.9.2014
% to do the simulation and get the wealth distribution.

% Preparing work
    clc
    clear
    close all

% Read data and normalize
     P = csvread('DataIndices.csv',1,1,[1,1,2530,2]);   % Daily Price data of DAX30, FTSE100 and SP500 from 12.31.2004 to 11.9.2014
 [r,c] = size(P);
     R = zeros(r-1,c);

for i = 2:r
    R(i-1,:) = log(P(i,:)./P(i-1,:));        % transfer the price to return data
end

% Generate the simulated data of 
       s = rng;                                     % set seeds
       n = 20000;                                % number of simulations
    mu = mean(R);           
sigma = cov(R);
  r_log = mvnrnd(mu,sigma,n);         % generate log returns
  r_dis = log(r_log+1);                       % discrete returns
 
 % Generate the wealth distribution
    W0 = 100;     % starting wealth
    fraction_interval = 0:1:100;
    W_T = zeros(n,length(fraction_interval));

for f = fraction_interval
    W_T(:,f+1) = W0*(1+[f/100,1-f/100]*r_dis');                     % Thorp Criterion to get the wealth distribution
end
   
    W_log = log(W_T/W0);
  W_Elog = mean(W_log);
    PW = max(W_Elog(mean(W_log <= -0.03) <= 0.01));     % finding the optimal returns based on constrains 
    num = find(W_Elog == PW)*0.01;
    opfrac = [num,1-num];                                                     % optimal fractions