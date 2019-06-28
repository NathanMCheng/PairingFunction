function [fit,resnorm] = ChiPeakFit(fun,fit,ydata,iterations)
%MMPFit Fits optical data using the MMP method

[fit,resnorm] = ...
    fmincon(@(x)norm(fun(x)-ydata),fit,[],[],[],[], ...
    [0,0,0],[inf,0.9067,0.4589],[],optimset('DiffMinChange',1e-4,'MaxFunEvals',iterations,'TolFun',1e-10));

end
