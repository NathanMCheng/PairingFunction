function [cond] = MemtoCond(mem1,mem2,w)
%ConvertMemoryFunctiontoConductivity 
%wp = 1; so do not calculate

opticalSelfE = (mem1+1i*mem2);
mu0 = 4*pi()*1e-7;
c = 2.99792e8;
kk = c*mu0/(2*pi());

cond = 1i/kk/(4*pi())./(w+opticalSelfE);


end

