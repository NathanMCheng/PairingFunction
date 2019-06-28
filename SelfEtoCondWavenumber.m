function [cond,w,integrand2] = SelfEtoCondWavenumber(selfE,w,T,impScat,wp)
%SelfEnergytoConductivity converts a given self energy and impurity
%scattering rate into a conductivity for corresponding frequency and temperature
%   Given a  self energy (selfE), impurity scattering rate (impScat) and
%   plasma frequency
%   for temperature (T) scalar and frequency (w) vector, integrates with respect to the
%   dielectric function (epsilon)
%   (nf(w+epsilon,T)-nf(epsilon,T))/(w-selfE(epsilon+w,T)+selfE*(epsilon,T)+1i*impScat)
%   where nf is the fermi dirac distribution and then multiplies by
%   wp^2/(4ipiw) where wp is the plasma frequency

%kb = 1.3806488e-23; %boltzmann constant
cond = ones(size(w)); %initiate conductivity
i = 0; %initate index of conductivity
Ktocminv = 0.695;
beta = (T/Ktocminv)^(-1);

% spacing = w(2)-w(1);
for wi = w
    i = i+1;
    integrand1 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
        +1).^(-1))/(wi+1i*impScat);
    
    k = 0; %initiate index corresponding to selfE(k) where k corresponds to epsilon
    integrand2 = zeros(1,length(w)); %initiate integrand for nonzero selfE epsilon>0
    integrand3 = zeros(1,length(w));
    for ep = -w
        k = k+1;
        
        if (k-i>0) 
            integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
                +1)^(-1))/(wi+selfE(k-i)-conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
        elseif (k==i)
            integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
                +1)^(-1))/(wi-conj(selfE(k))+1i*impScat);
        elseif (k-i<0) 
            integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
                +1)^(-1))/(wi-selfE(i-k)-conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e>wmax
        
        end
    end
    k = 0;
    for eps = w
        k = k+1;
        
        if (k+i<=length(selfE)) 
            integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
        +1)^(-1))/(wi-selfE(k+i)+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
        
        else
            integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
        +1)^(-1))/(wi+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e>wmax
        end
    end

    integrand4 = @(epsilon) ((exp(beta*(wi+epsilo    k = 0; %initiate index corresponding to selfE(k) where k corresponds to epsilonn))+1).^(-1)-(exp(beta*epsilon) ...
        +1).^(-1))/(wi-selfE(i)-imag(selfE(1))+1i*impScat);
    integrand5 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
        +1).^(-1))/(wi+selfE(i)-conj(selfE(1))+1i*impScat);
%     integrand3 = zeros(1,i-1); %initiate integrand for nonzero selfE epsilon<0
%     
%     k=0;
%     for ep = -w(1:i-1)
%         k=k+1;
%         integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%         +1)^(-1))/(wi+selfE(-k+i)+1i*impScat);
%     end
    
    cond(i) = wp^2/(4i*pi*wi)*( ... %(trapz(w,integrand2) + ... %Calculation of cond as a function of w,T for piecewise integral (selfE nonzero epsilon>0)
        trapz(w,integrand3) + ... %); % + ... %selfE nonzero epsilon<0
        integral(integrand1,-inf,-max(w)) + ... %selfE zero
        integral(integrand1,max(w),inf) + ...
        integral(integrand1,-min(w),min(w)) + ...
        integral(integrand4,0,min(w)) + ...
        integral(integrand5,-min(w),0)); %selfE zero
    

end
end


