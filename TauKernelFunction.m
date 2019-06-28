function [inverseTau] = TauKernelFunction(w,T,Pi)
%KERNELFUNCTION Calculates the kernel function for given w (eV) and T

KpereV = 11604.505;

if iscolumn(w)
    w = w';
end

% hcotangent1 = coth(w/(2*T/KpereV));
denominator = 2*pi()*T/KpereV;
L = zeros(length(w),length(w));
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
    L(i,:) = 1./w.*(2*w.*coth(wprime/denominator) - (w+wprime).*coth((w+wprime*ones(size(w)))/denominator) ... 
        + (w-wprime).*coth((w-wprime*ones(size(w)))/denominator));    
    
end
for k = 1:length(w)
    L(k,k) = pi()/w(k)*(2*w(k)*coth(w(k)/denominator)-2*w(k)*coth(2*w(k)/denominator));
end
% L(1,1) = L(2,1);
% for k = 2:length(w)
%     L(k,k) = L(k-1,k);
% end

integrand = zeros(size(L));
for k = 1:length(w)
    integrand(:,k) = Pi(w)'.*L(:,k);
end
inverseTau = trapz(w,integrand,1);

end
