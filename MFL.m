function [Pi] = MFL(A,T,wc,delta)
%MarginalFermiLiquid Computes the Bosonic spectral function via the
%Marginal Fermi Liquid theory
%   Detailed explanation goes here
KpereV = 11604.505;
Pi = @(w) A*1e7*tanh(w/(2*T*KpereV))*1./(1+exp((w-wc)/delta));

end

