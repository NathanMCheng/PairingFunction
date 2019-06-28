function [tau] = FittoTauX(x,w,kernel,T,dstart,dend)
%PitoDielectric Converts a bosonic spectral density into a tau function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec, dielec to
%   scattering rate
impScat = 0;
wp = 1;

Pi = Piecewise(x);
selfE = KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
dielec = CondtoDielec(cond,w);
tau = DielectoScattering(w,dielec,1);
tau = tau(dstart:dend);
end