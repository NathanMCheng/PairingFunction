function [fit,resnorm] = MFLFit(fun,fit,ydata,iterations)
%MMPFit Fits optical data using the MMP method

[fit,resnorm] = ...
    fmincon(@(x)norm(fun(x)-ydata),fit,[],[],[],[], ...
    [0,0,0,1,0,1],[inf,inf,inf,5,0.2,5],[],optimset('DiffMinChange',1e-1,'MaxFunEvals',iterations,'Algorithm','interior-point'));

end
