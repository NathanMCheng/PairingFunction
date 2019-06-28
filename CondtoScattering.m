function [inverseTau] = CondtoScattering(cond)
%CondToDielec turns conducitivy data into dielectric data
mu0 = 4*pi()*1e-7;
c = 2.99792e8;
kk = c*mu0/(2*pi());
cond = kk*cond;
inverseTau = 1/(4*pi())*real(cond)./(real(cond).^2+imag(cond).^2);

end

