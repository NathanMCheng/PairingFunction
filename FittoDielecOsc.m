function epsilon = FittoDielecOsc(x,w,kernel,T,dstart,dend,w0,coupling)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = 1;

Pi = Piecewise(x);
selfE = KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCond(selfE,w,T,impScat,wp);
epsilon = CondtoDielecNegative(cond,w) + LorentzOscill(w,w0,wp,coupling);
epsilon = epsilon(dstart:dend,:);
end
