% Title: A Varying-Coefficient Expectile Model
% Date: Dec 15, 2015
% Author: Dingshi Tian, Kirill Efimov

% clearing work&preparing
clc
clear
close all

% parameter setting
 M     = 500;                % number of simulations
 C     = 1;                  % choose the different DGP
 ng    = 5;                  % number of initial grid point
 nreg  = 3;                  % number of regressors(including 1!!!)
 gamma = 0.25;               % probability level
 TT    = [200 400 800];      % time intervals 
 RASE  = zeros(M, nreg); 
 RASE0 = zeros(2*length(TT), nreg);

% monte carlo simulations for M times
for n = 1:length(TT)
    T = TT(n);
    parfor j = 1:M
        [y, X] = mydgp(T, C);                    % data generating process
        if C == 1   
            u = X(:, 2);                          % smooth variable
        elseif C == 2
            u = X(:, 3);
        else 
            u = X(:, 4);
            X = X(:, 1:3);
        end
        grid = linspace(0.8*min(y), 0.8*max(y), 100)';
        b    = iniestor(y, X, u, grid, T);  
        beta = zeros(2*nreg-1, length(grid));                     
   
%  Estimate the functional coeffient of all the grid point
       % some lines for plotting
       % bp = zeros(1,length(u0));
       % bp0 = bp;  
       for i = 1:length(grid)
           beta(:, i) = myest(y, X, b(:, i), u, grid(i), T, gamma);  % use IWLLS to calculate the parameters for initial point
       end
                
%%  Calculate the RASE(root average squared error)
       if C == 1
           beta0 = [zeros(length(grid), 1), (0.138 + (0.316 + 0.982*grid).*exp(-3.89*grid.^2)),...
                   ...-0.437 - (0.659 + 0.126*grid).*exp(-3.89*grid.^2)]';  % real beta   
       elseif C == 2
           beta0 = [zeros(T, 1), 0.4*(U <= 1) - 0.8*(U > 1), -0.6*(U <= 1) + 0.8*(U > 1)]';
       elseif C == 3
           beta0 = [zeros(T, 1), sin(sqrt(2)*pi*U), cos(sqrt(2)*pi*U)]';
       else
           disp('No such DGP') 
       end
       RASE(j, :) = (1/T*sum((beta(1:nreg, :) - beta0).^2, 2)).^0.5;            % Calculate the root average squared error
    end
   
%% Store the parameters of each simulation
    RASE0(2*n - 1, :) = median(RASE);
    RASE0(2*n, :)     = std(RASE); 
end

% plot(u0,bp0,'k-',u0,bp,'k.')
