function Reflec = ChiFittoReflecX(x,chi,w,kernel,T,fitw)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
wp = x(2);
einf = x(3);

% Pi = griddedInterpolant(w,imag(chi(w)));
Pi = @(x) chi(x);
selfE = x(1)/250*KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
cond = griddedInterpolant(w,cond/60);
cond = cond(w);
dielec = CondtoDielec(transpose(cond),w,einf);
Reflec = DielectoRef(real(dielec),imag(dielec))';
Refi = griddedInterpolant(w,Reflec);
Reflec = Refi(fitw);
% Reflec = Reflec(1:375);
end
