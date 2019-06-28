function [L] = KernelFunction(w,T)
%KERNELFUNCTION Calculates the kernel function for given w (eV) and T

KpereV = 11604.505;

if iscolumn(w)
    w = w';
end

hcotangent = -1i*pi()*coth(w./(2*T/KpereV));
% hcotangent = -2*pi()*(1./(exp(w/T/KpereV)-1)+1/2);

% hcotangent = log(sinh(w./(2*T/KpereV)))+log(2);
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


for wprime = w;
    i=i+1;
    L(i,:) = hcotangent + double(psi(sym(1/2+1i*(w-wprime*ones(size(w)))/denominator)) - ... 
        psi(sym(1/2-1i*(w+wprime*ones(size(w)))/denominator)));
    
    
end

% for wprime = w;
%     i=i+1;
%     L(i,:) = hcotangent + mypsi(1/2+1i*(w-wprime)/denominator) + ... 
%         mypsi(1/2-1i*(w+wprime)/denominator);
%     
%     
% end
% 
% L = L*-1i*2*pi()*T/KpereV;


% L = cotmatrix + double(psi(sym(psi1))) - double(psi(sym(psi2)));


end

