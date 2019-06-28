function [kernel] = KernelPi(T)
%PiKernelFunction Calculates the KernelFunction relating bosonic
%excitations to optical conductivity and ellipsometry for a given
%Temperature T in kelvin.

m = 4001;
n = 500;

%Initiate
kernel = zeros(m,n);

%Unit Conversions
KpereV = 11604.505;

dw = 1e-3;

%X(w)
x0 = -2:dw:2; %w
x1 = -2.5:dw:2-dw; % w- w'
x2 = -2+dw:dw:2.5; %w + w'

%Terms
hcotangent = -1i*pi()*coth(x0./(2*T/KpereV));
denominator = 2*pi()*T/KpereV;

tic
psi1 = double(psi(sym(1/2+1i*x1/denominator)));
psi2 = double(psi(sym(1/2-1i*x2/denominator)));
toc

tic
for i = 1:m % w (-2,2)
    for j = 1:n %w' (0,0.5)
        kernel(i,j) = hcotangent(i)+psi1(500+i-j)-psi2(i+j-1);
    end
end
toc        
end

