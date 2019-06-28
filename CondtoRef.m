function [R,w] = CondtoRef(sigma,w,einf)
%   Converts conductivity (sigma) to reflectance (R)
%   Converts a complex conductivity to a reflectivity coefficient and
%   calculates the reflectance
mu0 = 4*pi()*1e-7;
c = 2.99792e8;
kk = c*mu0/(2*pi());
r=zeros(size(sigma)); %initialize reflectivity coefficient (r)
R=zeros(size(sigma)); %initialize reflectance
for i = 1:length(sigma)
    epsilon = kk*1i*4*pi()*sigma(i)/w(i)+einf; %dielectric function
    r(i) = ((epsilon^(1/2)-1)/(epsilon^(1/2)+1)); %calculation of reflectivity coefficient
    R(i) = r(i)*conj(r(i)); %calculation of reflectance

end

