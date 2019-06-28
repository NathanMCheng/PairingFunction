function [kernel] = KernelFunctionNegative(w,T)
%KERNELFUNCTION Calculates the kernel function for given w (eV) and T

KpereV = 11604.505;

if iscolumn(w)
    w = w';
end

cotangentw = -1i*pi()*coth(w/(2*T/KpereV));
cotangentww = -1i*pi()*coth(-w/(2*T/KpereV));
denominator = 2*pi()*T/KpereV;
L = zeros(length(w),length(w));
M = zeros(length(w),length(w));
% psi1 = zeros(length(w),length(w));
% psi2 = zeros(length(w),length(w));
i = 0; %initiate index
% cotmatrix = zeros(length(w),length(w));
% for wprime = w
%     i = i+1;
%     cotmatrix(i,:) = cotangent;
%     psi1(i,:) = (1/2+1i*(w-wprime*ones(size(w)))/denominator);
%     psi2(i,:) = (1/2-1i*(w+wprime*ones(size(w)))/denominator);
% end

for wprime = w
    i=i+1;
    L(i,:) = cotangentw + double(psi(sym(1/2+1i*(w-wprime*ones(size(w)))/denominator)) - ... 
        psi(sym(1/2-1i*(w+wprime*ones(size(w)))/denominator)));
    M(i,:) = cotangentww + double(psi(sym(1/2+1i*(-w-wprime*ones(size(w)))/denominator)) - ... 
        psi(sym(1/2-1i*(-w+wprime*ones(size(w)))/denominator)));
        
end
M = fliplr(M);
% L = cotmatrix + double(psi(sym(psi1))) - double(psi(sym(psi2)));
kernel = [M,L];

end

