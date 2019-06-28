function [fit,resnorm] = MMPFitTrust(fun,fit,xdata,ydata,iterations)
%MMPFit Fits optical data using the MMP method

[fit,resnorm] = ...
    lsqcurvefit(fun,fit,xdata,ydata,[0,0.01],[inf,0.21],optimset(... 
    'MaxFunEvals',iterations,'DiffMinChange',1e-1));

end