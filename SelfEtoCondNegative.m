function [cond,int] = SelfEtoCondNegative(selfE,fitww,T,impScat,wp)
%SelfEnergytoConductivity converts a given self energy and impurity
%scattering rate into a conductivity for corresponding frequency and temperature
%   Given a  self energy (selfE), impurity scattering rate (impScat) and
%   plasma frequency
%   for temperature (T) scalar and frequency (w) vector, integrates with respect to the
%   dielectric function (epsilon)
%   (nf(w+epsilon,T)-nf(epsilon,T))/(w-selfE(epsilon+w,T)+selfE*(epsilon,T)+1i*impScat)
%   where nf is the fermi dirac distribution and then multiplies by
%   wp^2/(4ipiw) where wp is the plasma frequency
if iscolumn(fitww)
    fitww = fitww';
end
% mu0 = 4*pi()*1e-7;
% c = 2.99792e8;
% kk = c*mu0/(2*pi());
%kb = 1.3806488e-23; %boltzmann constant
cond = ones(size(fitww)); %initiate conductivity
i = 0; %initate index of conductivity
KpereV = 11604.505;
kbt = T/KpereV;
% beta = (T/KpereV)^(-1);
% ww = -ones(size(w))*max(w)-w;
doubleselfE = [-conj(fliplr(transpose(selfE))),transpose(selfE)];
dE = [-fliplr(fitww),fitww];
% interpselfE = griddedInterpolant(fitww,selfE,'nearest');
interpselfE = griddedInterpolant(dE,doubleselfE,'linear');
% spacing = w(2)-w(1);
int = zeros(length(fitww),200);
for wi = fitww
    i = i+1;
    epsmin = -fitww(i)-21*kbt;
    epsmax = 21*kbt;
    de = (epsmax-epsmin)/200;
    w = epsmin:de:epsmax-de;
    ew = (wi+w)./kbt;
    e = w./kbt;
%     nfew = 1./(exp(ew)+1);
%     nfe = 1./(exp(e)+1);
nfew = zeros(length(w));
nfe = zeros(length(w));
    for m = 1:length(w)
        if (abs(ew(m))<=60)
            nfew(m) = 1./(exp(ew(m))+1);
        elseif ew(m)>0
            nfew(m) = 0;
        else nfew(m) = 1;
        end
        if (abs(e(m))<=60)
            nfe(m) = 1./(exp(e(m))+1);
        elseif e(m)>0
            nfe(m) = 0;
        else nfe(m) = 1;
        end
    end
%     integrand1 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
%         +1).^(-1))/(wi-interpselfE(wi+epsilon)+conj(interpselfE(epsilon)) +1i*impScat);
%     
%     k = 0; %initiate index corresponding to selfE(k) where k corresponds to epsilon
%     integrand2 = zeros(1,length(w)); %initiate integrand for nonzero selfE epsilon>0
%     for ep = -w
%         k = k+1;
%         
%         if (ep+wi<0) 
%             integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%                 +1)^(-1))/(wi+interpselfE(-(wi+ep))-conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
%         elseif (ep+wi>=0) 
%             integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%                 +1)^(-1))/(wi-interpselfE(wi+ep)-conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e>wmax        
%         end
%     end
%     integrand6 = zeros(1,length(w));
%     k = 0;
%     for epsi = ww
%         k = k+1;
%         integrand6(k) = ((exp(beta*(wi+epsi))+1).^(-1)-(exp(beta*epsi) ...
%         +1).^(-1))/(wi+interpselfE(-(wi+epsi))+1i*impScat);
%     end 
    
    
    integrand3 = zeros(1,length(w));
    k = 0;
    selfE1 = interpselfE(w);
    selfE2 = interpselfE(wi+w);
%     selfE1 = (real(selfE1)*kbt+1i*imag(selfE1);
%     selfE2 = (real(selfE2)*kbt+1i*imag(selfE2);
    for eps = w
        k = k+1;
        integrand3(k) = (nfew(k)-nfe(k))./(wi-selfE2(k)+conj(selfE1(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
    end
%         integrand4 = zeros(1,500);
%     k = 0;
%     selfE1 = interpselfE(w);
%     selfE2 = interpselfE(wi+w);
%     for eps = w(501:1000)
%         k = k+1;
%         integrand4(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
%             +1)^(-1))/(wi-selfE2(k)+conj(selfE1(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
%     end


%     integrand4 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
%         +1).^(-1))/(wi-selfE(i)-imag(selfE(1))+1i*impScat);
%     integrand5 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
%         +1).^(-1))/(wi+selfE(i)+imag(selfE(1))+1i*impScat);
%     integrand3 = zeros(1,i-1); %initiate integrand for nonzero selfE epsilon<0
%     
%     k=0;
%     for ep = -w(1:i-1)
%         k=k+1;
%         integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%         +1)^(-1))/(wi+selfE(-k+i)+1i*impScat);
% %     end
%     if mod(i,2)==0
%         plot(w,imag(integrand3),'r');
%     else plot(w,imag(integrand3),'b');
%     end
% 
%     hold on; 
%     int(i,:) = integrand3;
    
    cond(i) = wp^2/(4i*pi()*wi)*(... %trapz(w,integrand2) + ... %Calculation of cond as a function of w,T for piecewise integral (selfE nonzero epsilon>0)
        + ... trapz(w(501:1000),integrand4) + ...
        de*trapz(integrand3)) ; %+ ... %); % + ... %selfE nonzero epsilon<0
        %...%integral(integrand1,-inf,-2*max(w)) + ... %selfE zero
        %... integral(integrand1,max(w),inf,'ArrayValued',true) + ... % + ...
        %... %... %integral(integrand1,-min(w),min(w)) +...
        %integral(integrand1,0,inf,'ArrayValued',true)); %+ ...
        %integral(integrand5,-min(w),0)); %selfE zero
    

end
end
