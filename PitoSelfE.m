function [selfE,L,integrand] = PitoSelfE(Pi,w,T)
%PItoSelfEnergy converts a histrogramic Pi function to a self energy for
%temperature (T) and frequency (w) in eV
%   Integrates Pi(w')L(w,w',T) over frequency to give self energy
%   selfE(w,T) where L(w,w',T) = -i*pi*coth(w/2T) + digamma(1/2+i(w-w')/2piT) -
%   digamma(1/2-i(w+w')/2piT) and digamma is the digamma function.
%   w is a vector, T is a scalar

%selfE = ones(1,11); %initiate self energy


% hbar = 1.05457173e-34; %planck's constant
% eV = 1.60217657e-19; %one eV in J

% w = w/hbar*eV;
%kb = 1.3806488e-23; %boltzmann constant
KpereV = 11604.505;
beta = (T/KpereV)^(-1);

L = zeros(length(w),length(w));
integrand = zeros(length(w),length(w));
i = 0; %initiate index
for wprime = w;
    i=i+1;
    L(i,:) = -2i*pi()*((exp(beta*w)-1).^(-1)+1/2)+ double(psi(sym(1/2+1i*(w-wprime*ones(size(w)))/(2*pi()*T/KpereV))) - ... 
        psi(sym(1/2-1i*(w+wprime*ones(size(w)))/(2*pi()*T/KpereV))));
    
    
end

for k = 1:i
    integrand(:,k) = Pi(w)'.*L(:,k);
end


selfE = trapz(w,integrand);

end





% for wi = w
%     i = i+1;
%     
%     L = @(wprime) -1i*pi()*coth(wi/2*T) + double(psi(sym(1/2+1i*(wi-wprime)/(2*pi()*T))) - psi(sym(1/2-1i*(wi+wprime)/(2*pi()*T)))); % defines L
%     integrand = @(x) L(x).*Pi(x); %defines the integral
% 
%     selfE(i) = integral(integrand,0,inf); %integrates over wprime to solve for self energy
% 
% end

%-2i*pi()*(exp(w*(kb*T)^(-1)-1).^(-1)+1/2) science
%-i*pi()*coth(w/(2*T)) van heumen