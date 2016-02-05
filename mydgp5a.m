function[a,b]    =  mydgp5a(T)
 
uu = zeros(T+100,1);
deltau = zeros(T+100,1);
Y = zeros(T+102,1);
    
for t = 1:(T+100)
        
    uu(t) = unifrnd(0,3);
    deltau(t) = 3*exp(-4*(uu(t)-1)^2)+2*exp(-5*(uu(t)-2)^2);
    a1 = sin(sqrt(2)*pi*uu(t));
    a2 = cos(sqrt(2)*pi*uu(t));
    Y(t+2) = a1*Y(t+1) + a2*Y(t) +deltau(t)*normrnd(0,1);
        
end
    
a = Y(103:T+102);                               % Y(t):delete the first 100 values 
b = [ones(T,1),Y(102:T+101),Y(101:T+100),uu(101:T+100)];      % X(t) = [1, Y(t-1), Y(t-2),U]