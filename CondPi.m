function [cond] = CondPi(selfE,T,wp)
%ConductivityFunction Calculates the conductivity from the self energy as
%calculated from the bosonic spectrum. T is temperature in Kelvin. wp is
%the plasma frequency. selfE is the self energy.

%Unit Conversions
KpereV = 11604.505;
kbt = T/KpereV;
kbt21 = 21*kbt;
dw = 1e-3;
wp2 = wp^2;

%Initialize
cond = zeros(1000,1);

%X(w)
x1 = dw:dw:1; %w
selfEw = -2:dw:2; %w

selfEinterp = griddedInterpolant(selfEw,selfE);

for i = 1:1000 %w
    
    %Terms
    de = (kbt21-(-kbt21-x1(i)))/200;
    
    x0 = -kbt21-x1(i):de:kbt21; %e
    x2 = -kbt21:de:kbt21+x1(i); %e +w
    
    nfe = (exp(x0/kbt)+1).^(-1); %nf(e)
    nfew = (exp(x2/kbt)+1).^(-1); %(nf(e+w))
    
    selfEe = selfEinterp(x0);
    selfEew = selfEinterp(x2);
    
    
%     cond(i) = wp2/(4i*pi()*x1(i))*...
%         de*trapz((nfew-nfe)./ ...
%         (x1(i)-selfEew+conj(selfEe)));
plot(imag((nfew-nfe)./ ...
        (x1(i)-selfEew+conj(selfEe))));
    hold on
    
end

