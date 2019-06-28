function epsilon = MMPFittoDielecX(x,w,kernel,T,dstart,dend)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = 1;

Pi = MMP(x(1),x(2),x(3));
selfE = KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
epsilon = CondtoDielec(cond,w);
epsilon = epsilon(dstart:dend,:);
end
