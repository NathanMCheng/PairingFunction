function [fit,resnorm] = MMPFitLevenberg(fun,fit,xdata,ydata,iterations)
%MMPFit Fits optical data using the MMP method

[fit,resnorm] = ...
    lsqcurvefit(fun,fit,xdata,ydata,[],[],optimset('Algorithm','levenberg-marquardt', ... 
    'MaxFunEvals',iterations));

end