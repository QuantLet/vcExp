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
 
function[a, b] = mydgp(T, C)           % sub-functions for data generating

% matrix setting 
Y    = zeros(T + 102, 1);
Y(2) = 1;
          
% generate initial value of Y  
if C == 1                               % Example 1
    for t = 1 : (T + 100)
        uu       = Y(t + 1);
        a1       = 0.138 + (0.316 + 0.982*uu)*exp(-3.89*uu^2);
        a2       = -0.437 - (0.659 + 0.126*uu)*exp(-3.89*uu^2);
        eps      = normrnd(0, 0.2);
        Y(t + 2) = Y(t + 1)*a1 + Y(t)*a2 + 2*a1*Y(t + 1)*eps*(uu >= 0);
    end
    a = Y(103 : T + 102);                               % Y(t):delete the first 100 values
    b = [ones(T, 1), Y(102 : T + 101), Y(101 : T + 100)];      % X(t) = [1, Y(t-1), Y(t-2)]
elseif C == 2                            % Example 2
    for t = 1 : (T + 100)
        uu       = Y(t);
        a1       = 0.4*(uu <= 1) - 0.8*(uu > 1);
        a2       = -0.6*(uu <= 1) - 0.2*(uu > 1);
        Y(t + 2) = Y(t + 1)*a1 + Y(t)*a2 + normrnd(0, 1);
    end
    a = Y(103 : T + 102);                               % Y(t):delete the first 100 values
    b = [ones(T, 1), Y(102 : T + 101), Y(101 : T + 100)];      % X(t) = [1, Y(t-1), Y(t-2)]
elseif C == 3                             % Example 3
    uu = zeros(T + 100, 1);
    for t = 1 : (T + 100)
        uu(t) = unifrnd(0, 3);
        a1    = sin(sqrt(2)*pi*uu(t));
        a2    = cos(sqrt(2)*pi*uu(t));
        r     = unifrnd(0, 1);
        if r >= 0 && r < 0.95
            eps = normrnd(0, 1/sqrt(0.95));
        else
            eps = normrnd(-5, 1/sqrt(0.05));
        end
        Y(t + 2) = a1*Y(t + 1) + a2*Y(t) + eps;
    end
    
    a = Y(103 : T + 102);                                             % Y(t):delete the first 100 values
    b = [ones(T, 1),Y(102 : T + 101), Y(101 : T + 100), uu(101 : T + 100)];      % X(t) = [1, Y(t-1), Y(t-2),U]
end

function[b]=myest(y, X, u, ug, T)              % Obtaining initial estimators

% Variable setting 
bd02   = 1.06*std(y)*T^(-0.2);
[~, c] = size(X);
b      = zeros(2*c - 1, length(ug));
S      = zeros(2*c - 1, 2*c - 1);
R      = zeros(2*c - 1, 1);
 
% obatain the initial value of estimators 
for i = 1:length(ug)
    z = [X(:, 1), X(:, 2), X(:, 3), X(:, 2).*(u - ug(i)), X(:, 3).*(u - ug(i))];
    for t = 1:T
        kernel01 = 1/sqrt(2*pi)*exp(-0.5*((u(t) - ug(i))/bd02)^2);
        s        = kernel01*z(t, :)'*z(t, :);
        S        = s + S;
        r        = kernel01*z(t, :)'*y(t);
        R        = r + R;
    end
    b(:, i) = S\R;
end

return(b)

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
        beta = zeros(2*nreg - 1, length(grid));                     
   
%  Estimate the functional coeffient of all the grid point
       % some lines for plotting
       % bp = zeros(1,length(u0));
       % bp0 = bp;  
       for i = 1:length(grid)
           beta(:, i) = myest(y, X, b(:, i), u, grid(i), T, gamma);           % use IWLLS to calculate the parameters for initial point
       end
                
%%  Calculate the RASE(root average squared error)
       if C == 1
           beta0 = [zeros(length(grid), 1), (0.138 + (0.316 + 0.982*grid).*exp(-3.89*grid.^2)),...
                   ...-0.437 - (0.659 + 0.126*grid).*exp(-3.89*grid.^2)]';     % real beta   
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
