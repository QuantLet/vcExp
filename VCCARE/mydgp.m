function[a, b] = mydgp(T, C)

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
