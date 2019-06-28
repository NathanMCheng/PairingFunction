function [selfE] = PitoSelfEsym(Pi)
%PitoSelfE converts Pi(w) to selfE(w) symbolically
%   symbolically integrates Pi(w')L(w,w',T) over frequency to give self energy
%   selfE(w,T) where L(w,w',T) = -i*pi*coth(w/2T) + digamma(1/2+i(w-w')/2piT) -
%   digamma(1/2-i(w+w')/2piT) and digamma is the digamma function.
%   w is a vector, T is a scalar

syms w T wprime

L(w,wprime,T) = -1i*pi()*coth(w/2*T) + psi(sym(1/2+1i*(w-wprime)/(2*pi()*T))) - ... 
    psi(sym(1/2-1i*(w+wprime)/(2*pi()*T))); % defines L

selfE(w,T) = int(L(w,wprime,T)*Pi(wprime),wprime,0,inf);

%selfE = matlabFunction(selfE);

end

