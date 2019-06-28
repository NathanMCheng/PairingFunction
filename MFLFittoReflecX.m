function dielec = MFLFittoReflecX(x,w,kernel,T)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
einf = x(6);
impScat = x(5);
wp = x(4);

Pi = MFL(x(1),T,x(2),x(3));
selfE = KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
cond = griddedInterpolant(w,cond/60);
cond = cond(w);
dielec = CondtoDielec(transpose(cond),w,einf);
% Reflec = DielectoRef(real(dielec),imag(dielec))';
% Reflec = Reflec(1:359);
end
