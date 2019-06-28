function epsilon = ChiCutoffFittoDielecX(x,chi,w,kernel,T,dstart,dend)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = 1;

Pi = griddedInterpolant(w,(w<=x(2)).*imag(chi(w)));
selfE = (x(1)/1000)*KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
epsilon = CondtoDielec(cond,w);
epsilon = epsilon(dstart:dend,:);
end