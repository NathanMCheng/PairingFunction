function Dielec = rpaChiFittoDielecX(x,rpaChi,Chi,w,kernel,T,dstart,dend)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = 0.5;

selfE = x(1)*KernelPitoSelfE(rpaChi,kernel,w)+x(2)*KernelPitoSelfE(Chi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
Dielec = CondtoDielec(cond,w);
Dielec = Dielec(dstart:dend);
end
