function Reflec = FittoReflecX(x,w,kernel,T)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
m = length(x);
einf = x(m);
impScat = x(m-1);
wp = x(m-2);

Pi = Piecewise(x(1:m-3));
selfE = KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
cond = griddedInterpolant(w,cond/60);
cond = cond(w);
dielec = CondtoDielec(transpose(cond),w,einf);
Reflec = DielectoRef(real(dielec),imag(dielec))';
Reflec = Reflec(1:359);
end
