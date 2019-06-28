function Reflec = rpaChiFittoReflecX(x,rpaChi,w,kernel,T)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = x(2);
einf = x(3);

Pi = griddedInterpolant(w,(rpaChi(w)));
selfE = x(1)*KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
cond = griddedInterpolant(w,cond/60);
cond = cond(w);
dielec = CondtoDielec(transpose(cond),w,einf);
Reflec = DielectoRef(real(dielec),imag(dielec))';
Reflec = Reflec(1:375);
end
