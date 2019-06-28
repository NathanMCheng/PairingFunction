function [epsi] = CondtoDielec(cond,w,einf)
%CondToDielec turns conducitivy data into dielectric data
if iscolumn(w)
    w = w';
end
mu0 = 4*pi()*1e-7;
c = 2.99792e8;
kk = c*mu0/(2*pi());
epsi = 4i*pi()*kk*cond./w+einf;
%epsi = [real(epsi)',imag(epsi)'];
epsi = transpose(epsi);
end

