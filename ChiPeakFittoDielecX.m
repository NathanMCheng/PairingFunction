function epsilon = ChiPeakFittoDielecX(x,chi,w,kernel,T,dstart,dend)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = 1;

Chi = griddedInterpolant(w,imag(chi(w)));
%LowEPeak = Piecewise([0.0138,0.0085,x(2),0.0190,x(3),0.0539]);
LowEPeak = Piecewise([0,0.01,0.04,0.04,0.4,0.09]);
%selfE = x(1)*KernelPitoSelfE(Chi,kernel,w) + KernelPitoSelfE(LowEPeak,kernel,w);
selfE = x*KernelPitoSelfE(Chi,kernel,w) + KernelPitoSelfE(LowEPeak,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
epsilon = CondtoDielec(cond,w);
epsilon = epsilon(dstart:dend,:);
end
