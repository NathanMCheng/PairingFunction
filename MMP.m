function [Pi] = MMP(coupling,gamma,cutoff)
%Millis-Monien-Pines Computes the Miliens-Monien-Pines model for the
%bosonic spectral function

Pi = @(w) (w<cutoff)*coupling.*w./(gamma^2+w.^2);

end

