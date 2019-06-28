function [fit,resnorm] = GenericFit(fun,fit,ydata,iterations)
%MMPFit Fits optical data using the MMP method

[fit,resnorm] = ...
    fmincon(@(x)norm(fun(x)-ydata),fit,[],[],[],[], ...
    [0,5e-4,0,0],[0.005,inf,inf],[],optimset('DiffMinChange',1e-4,'MaxFunEvals',iterations));

end
