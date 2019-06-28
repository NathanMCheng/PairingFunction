function [selfE,integrand] = GIKernelPitoSelfE(Pi,Kernel,w)
%PItoSelfEnergy converts a histrogramic Pi function to a self energy for
%temperature (T) and frequency (w) in eV
%   Integrates Pi(w')L(w,w',T) over frequency to give self energy
%   selfE(w,T) where L(w,w',T) = -i*pi*coth(w/2T) + digamma(1/2+i(w-w')/2piT) -
%   digamma(1/2-i(w+w')/2piT) and digamma is the digamma function.
%   w is a vector, T is a scalar

%selfE = ones(1,11); %initiate self energy
if iscolumn(w)
    w = w';
end

% hbar = 1.05457173e-34; %planck's constant
% eV = 1.60217657e-19; %one eV in J

% w = w/hbar*eV;
%kb = 1.3806488e-23; %boltzmann constant
%KpereV = 11604.505;


% L = zeros(length(w),length(w));
integrand = zeros(size(Kernel));
% i = 0; %initiate index
% for wprime = w;
%     i=i+1;
%     L(i,:) = -1i*pi()*coth(w/(2*T/KpereV))+ double(psi(sym(1/2+1i*(w-wprime*ones(size(w)))/(2*pi()*T/KpereV))) - ... 
%         psi(sym(1/2-1i*(w+wprime*ones(size(w)))/(2*pi()*T/KpereV))));
%     
%     
% end

for k = 1:length(Kernel)
    integrand(:,k) = Pi(w).*Kernel(:,k);
end


selfE = trapz(w,integrand,1);

end