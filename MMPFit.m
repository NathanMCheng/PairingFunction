function [fit,resnorm] = MMPFit(fun,fit,ydata,iterations)
%MMPFit Fits optical data using the MMP method

[fit,resnorm] = ...
    fmincon(@(x)norm(fun(x)-ydata),fit,[],[],[],[], ...
    [0,0.01,0.4,1,0,1],[inf,0.24,4.5,5,0.2,5],[],optimset('DiffMinChange',1e-1,'MaxFunEvals',iterations,'TolFun',1e-10,'Algorithm','interior-point'));

end
