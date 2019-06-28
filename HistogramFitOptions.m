function [fit,resnorm] = HistogramFitOptions(fun,fit,ydata,iterations)
%HistogramFit Fit's the dielectric function via histogramic bosonic
%spectral density

%   Designed to modify the lower bound continuously for the histogramic
%   representation

if length(fit) == 6
    [fit,resnorm] = ...
        fmincon(@(x)norm(fun(x)-ydata),fit,[0,-1,0,0,0,0;0,1,0,-1,0,0;0,0,0,1,0,-1],[0;0;0],[],[], ...
        [0,0,0,0,0,0],[],[],optimset('DiffMinChange',1e-2,'MaxFunEvals',iterations));
    
elseif length(fit) == 10    
    [fit,resnorm] = ...
        fmincon(@(x)norm(fun(x)-ydata),fit, ...
        [0,-1,0,0,0,0,0,0,0,0; ... 
         0,1,0,-1,0,0,0,0,0,0; ...
         0,0,0,1,0,-1,0,0,0,0; ...
         0,0,0,0,0,1,0,-1,0,0; ...
         0,0,0,0,0,0,0,1,0,-1],[0;0;0;0;0],[],[], ...
        [0,0,0,0,0,0,0,0,0,0,],[],[], ...
        optimset('DiffMinChange',1e-1,'MaxFunEvals',iterations,'TolX',1e-6));

elseif length(fit) == 12
    [fit,resnorm] = ...
        fmincon(@(x)norm(fun(x)-ydata),fit, ...
        [0,-1,0,0,0,0,0,0,0,0,0,0; ...
        0,1,0,-1,0,0,0,0,0,0,0,0; ...
        0,0,0,1,0,-1,0,0,0,0,0,0; ...
        0,0,0,0,0,1,0,-1,0,0,0,0; ...
        0,0,0,0,0,0,0,1,0,-1,0,0;
        0,0,0,0,0,0,0,0,0,1,0,-1],[0;0;0;0;0;0],[],[], ...
        [0,0,0,0,0,0,0,0,0,0,0,0],[],[], ...
        optimset('DiffMinChange',1e-1,'MaxFunEvals',iterations,'TolX',1e-6));
end

