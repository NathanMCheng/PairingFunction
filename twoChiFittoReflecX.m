function Reflec = twoChiFittoReflecX(x,rpaChi,Chi,w,kernel,T)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = x(3);
einf = x(4);

selfE = x(1)*KernelPitoSelfE(rpaChi,kernel,w)+x(2)*KernelPitoSelfE(Chi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
cond = griddedInterpolant(w,cond/60);
cond = cond(w);
dielec = CondtoDielec(transpose(cond),w,einf);
Reflec = DielectoRef(real(dielec),imag(dielec))';
Reflec = Reflec(1:375);
end
