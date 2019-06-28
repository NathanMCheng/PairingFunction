function [wChi] = IntegrateMomentum(wkqChi)
%IntegrateMomentum Integrates out the momentum dependence of a Chi/Pi
%function
% 

q = -pi():0.05:pi(); %q momentum
k = -pi():0.05:pi(); %k momentum

wChi = trapz(k,trapz(q,wkqChi,3),2);



end

