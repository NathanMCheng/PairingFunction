function [inverseTau] = DielectoScattering(frequency,epsilon,einfinity)
%DielecricFunctiontoScattering Converts dielectric function to scattering
%rate

inverseTau = imag(epsilon)./(frequency.*((einfinity- ... 
    real(epsilon)).^2+imag(epsilon).^2));


end

