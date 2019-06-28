function [cond] = SelfEtoCond(selfE,w,T,impScat,wp)
%SelfEnergytoConductivity converts a given self energy and impurity
%scattering rate into a conductivity for corresponding frequency and temperature
%   Given a  self energy (selfE), impurity scattering rate (impScat) and
%   plasma frequency
%   for temperature (T) scalar and frequency (w) vector, integrates with respect to the
%   dielectric function (epsilon)
%   (nf(w+epsilon,T)-nf(epsilon,T))/(w-selfE(epsilon+w,T)+selfE*(epsilon,T)+1i*impScat)
%   where nf is the fermi dirac distribution and then multiplies by
%   wp^2/(4ipiw) where wp is the plasma frequency

if iscolumn(w)
    w = w';
end
cond = ones(size(w)); %initiate conductivity
i = 0; %initate index of conductivity
KpereV = 11604.505;
beta = (T/KpereV)^(-1);
interpselfE = griddedInterpolant(w,selfE,'linear');
wp4ipi = wp^2/4i/pi();
for wi = w
    i = i+1;
    integrand3 = zeros(1,length(w));
    k = 0;
    selfE2 = interpselfE(wi+w);
    for eps = w
        k = k+1;
        
        if (eps+wi<=max(w)) 
            integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
        +1)^(-1))/(wi-selfE2(k)+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
        
        else
            integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
        +1)^(-1))/(wi+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e>wmax
        end
    end    
    cond(i) = wp4ipi/wi*(trapz(w,integrand3)) ;%selfE nonzero epsilon<0


end
end

















% %kb = 1.3806488e-23; %boltzmann constant
% cond = ones(size(w)); %initiate conductivity
% i = 0; %initate index of conductivity
% KpereV = 11604.505;
% beta = (T/KpereV)^(-1);
% mu0 = 4*pi()*1e-7;
% c = 2.99792e8;
% % kk = c*mu0/(2*pi());
% % spacing = w(2)-w(1);
% for wi = w
%     i = i+1;
%     integrand1 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
%         +1).^(-1))/(wi+1i*impScat);
%     
%     k = 0; %initiate index corresponding to selfE(k) where k corresponds to epsilon
%     integrand2 = zeros(size(w)); %initiate integrand for nonzero selfE epsilon>0
%     
%     for ep = w
%         k = k+1;
%         
%         if (k+i<=numel(selfE)) 
%             integrand2(k) = ((exp(beta*(wi+ep))+1).^(-1)-(exp(beta*ep) ...
%         +1).^(-1))/(wi-selfE(k+i)+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
%         
%         else 
%             integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%         +1)^(-1))/(wi+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e>wmax
%         end
%     end
%     
% %     integrand5 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
% %         +1).^(-1))/(wi+selfE(i)-conj(selfE(1))+1i*impScat);
% % 
% %     integrand3 = zeros(1,i-1); %initiate integrand for nonzero selfE epsilon<0
% %     
% %     k=0;
% %     for ep = -w(1:i-1)
% %         k=k+1;
% %         if (i-k>0)
% %             integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
% %                 +1)^(-1))/(wi-selfE(-k+i)-conj(selfE(k))+1i*impScat);
% %         elseif (i-k<0)
% %             integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
% %                 +1)^(-1))/(wi+selfE(-k+i)-conj(selfE(k))+1i*impScat);
% %         else
% %             integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
% %                 +1)^(-1))/(wi-conj(selfE(k))+1i*impScat);
% %         end
% %     end
% 
%     integrand4 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
%         +1).^(-1))/(wi-selfE(i)+conj(selfE(1))+1i*impScat);
%     
%     cond(i) = wp^2/(4i*pi*wi)*(trapz(w,integrand2)); % + ... %Calculation of cond as a function of w,T for piecewise integral (selfE nonzero epsilon>0)
%         %... %trapz(w(1:k),integrand3,2) + ... %selfE nonzero epsilon<0
%         %... %integral(integrand1,-inf,-w(i)) + ... %selfE zero
%         %integral(integrand1,max(w),inf) + ...
%         %integral(integrand4,0,min(w))); % + ... %selfE zero
%         %integral(integrand5,-min(w),0)); %selfE zero
%     
% 
% end
% end
% 
% % note to self: consider griddedInterpolant selfE and do the whole
% % integral.

