function [cond] = DielectoCond(w,epsilon,einf)
%DielectoConductivity Converts Dielectric Function data to Conductivity
mu0 = 4*pi()*1e-7;
c = 2.99792e8;
kk = c*mu0/(2*pi());
cond = 1i*w/kk/(4*pi()).*(einf-epsilon);


end

