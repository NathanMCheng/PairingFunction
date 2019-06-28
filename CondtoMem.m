function mem = CondtoMem(cond,w,wp,einf)
%ConductivityToMemoryFunction 
% wp = 1; do not calculate with it
if iscolumn(w)
    w = w';
end
if iscolumn(cond)
    cond = transpose(cond);
end
mu0 = 4*pi()*1e-7;
c = 2.99792e8;
kk = c*mu0/(2*pi());
mem = w.*((1i*wp^2./(w*kk*4*pi()))./(cond)-w.^2.*einf)-w;
% mem = -real(mem)+1i*imag(mem);
% mem1 = real(mem);
% mem2 = imag(mem);
end

