function [L] = Kernel(w,T)
%KERNELFUNCTION Calculates the kernel function for given w (eV) and T

KpereV = 11604.505;

if iscolumn(w)
    w = w';
end

bosedistribution = -2i*pi()*(exp(w/(T/KpereV)-1).^(-1)+1/2);
denominator = 2*pi()*T/KpereV;
L = zeros(length(w),length(w));
% psi1 = zeros(size(L));
% psi2 = zeros(size(L));
i = 0; %initiate index
% cotmatrix = zeros(size(L));
% for wprime = w
%     i = i+1;
%     cotmatrix(i,:) = cotangent;
%     psi1(i,:) = (1/2+1i*(w-wprime*ones(size(w)))/denominator);
%     psi2(i,:) = (1/2+1i*(w+wprime*ones(size(w)))/denominator);
% end

for wprime = w
    i=i+1;
    L(i,:) = bosedistribution + double(psi(sym(1/2+1i*(w-wprime*ones(size(w)))/denominator)) - ... 
        psi(sym(1/2-1i*(w+wprime*ones(size(w)))/denominator)));
    
    
end

% L = cotmatrix + double(psi(sym(psi1))) - double(psi(sym(psi2)));


end

