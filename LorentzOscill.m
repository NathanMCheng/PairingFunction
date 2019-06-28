function [epsilonL] = LorentzOscill(w,w0,wp,coupling,einf)
%LORENTZOSCILLATORS Calcualtes the lorentz dielectric function
%   calculation of the lorentz dielectric function takes into consideration
%   varying numbers of lorentz oscillators

epsilonL = einf*ones(size(w));

for i = 1:length(w0)
    epsilonL = epsilonL + wp(i)^2./(w0(i)^2*ones(size(w))-w.^2-1i*coupling(i)*w);
end


end

