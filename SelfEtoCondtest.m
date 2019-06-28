function [cond,dE] = SelfEtoCondtest(selfE,w,T,impScat,wp)
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
% mu0 = 4*pi()*1e-7;
% c = 2.99792e8;
% kk = c*mu0/(2*pi());
%kb = 1.3806488e-23; %boltzmann constant
cond = zeros(size(w)); %initiate conductivity
i = 0; %initate index of conductivity
KpereV = 11604.505;
beta = (T/KpereV)^(-1);
%ww = -ones(size(w))*max(w)-w;
doubleselfE = [-fliplr(selfE),selfE];
dE = [-fliplr(w),w];
positiveinterp = griddedInterpolant(w,selfE,'linear');
interpselfE = griddedInterpolant(dE,doubleselfE,'nearest');
negw = [-max(w)-fliplr(w),-fliplr(w)];
% spacing = w(2)-w(1);
for wi = w
    i = i+1;
    integrand1 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
        +1).^(-1))/(wi+1i*impScat);
    
    k = 0; %initiate index corresponding to selfE(k) where k corresponds to epsilon
    integrand2 = zeros(1,length(w)); %initiate integrand for nonzero selfE epsilon<0
    for ep = w
        k = k+1;
        
        integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
            +1)^(-1))/(wi-interpselfE(wi+ep)+conj(interpselfE(ep))+1i*impScat);
        
    end
    
    integrand3 = zeros(1,length(negw));
    
    for ep = negw
        integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
            +1)^(-1))/(wi-interpselfE(wi+ep)+conj(interpselfE(ep))+1i*impScat);
    end
    plot(negw,imag(integrand3))
    hold on
%     integrand6 = zeros(1,length(w));
%     k = 0;
%     for epsi = ww
%         k = k+1;
%         integrand6(k) = ((exp(beta*(wi+epsi))+1).^(-1)-(exp(beta*epsi) ...
%         +1).^(-1))/(wi+interpselfE(-(wi+epsi))+1i*impScat);
%     end 
% 
%     integrand3 = zeros(1,length(w));
%     k = 0;
%     for eps = w
%         k = k+1;
%         
%         if (eps+wi<=max(w)) 
%             integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
%         +1)^(-1))/(wi-interpselfE(wi+eps)+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
%         
%         else
%             integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
%         +1)^(-1))/(wi+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e>wmax
%         end
%     end

    integrand4 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
        +1).^(-1))/(wi-selfE(i)+conj(positiveinterp(epsilon))+1i*impScat);
    integrand5 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
        +1).^(-1))/(wi+selfE(i)-conj(positiveinterp(epsilon))+1i*impScat);
%     integrand3 = zeros(1,i-1); %initiate integrand for nonzero selfE epsilon<0
%     
%     k=0;
%     for ep = -w(1:i-1)
%         k=k+1;
%         integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%         +1)^(-1))/(wi+selfE(-k+i)+1i*impScat);
%     end
    
    cond(i) = wp^2/(4i*pi*wi)*(trapz(w,integrand2) + trapz(negw,integrand3) + ... %Calculation of cond as a function of w,T for piecewise integral (selfE nonzero epsilon>0)
        integral(integrand4,0,min(w),'ArrayValued',true) + ...
        integral(integrand5,-min(w),0,'ArrayValued',true) + ... %selfE zero    
        ... %... %trapz(ww,integrand6)); % + ...
        ... %trapz(w,integrand3)); % + ... %); % + ... %selfE nonzero epsilon<0
        integral(integrand1,-inf,-2*max(w),'ArrayValued',true) + ... %selfE zero
        integral(integrand1,max(w),inf,'ArrayValued',true)); % + ... % + ...
        %... %integral(integrand1,-min(w),min(w)) +...

    

end
end

% function [cond,dE] = SelfEtoCondtest(selfE,w,T,impScat,wp)
% %SelfEnergytoConductivity converts a given self energy and impurity
% %scattering rate into a conductivity for corresponding frequency and temperature
% %   Given a  self energy (selfE), impurity scattering rate (impScat) and
% %   plasma frequency
% %   for temperature (T) scalar and frequency (w) vector, integrates with respect to the
% %   dielectric function (epsilon)
% %   (nf(w+epsilon,T)-nf(epsilon,T))/(w-selfE(epsilon+w,T)+selfE*(epsilon,T)+1i*impScat)
% %   where nf is the fermi dirac distribution and then multiplies by
% %   wp^2/(4ipiw) where wp is the plasma frequency
% if iscolumn(w)
%     w = w';
% end
% mu0 = 4*pi()*1e-7;
% c = 2.99792e8;
% kk = c*mu0/(2*pi());
% %kb = 1.3806488e-23; %boltzmann constant
% cond = zeros(size(w)); %initiate conductivity
% i = 0; %initate index of conductivity
% KpereV = 11604.505;
% beta = (T/KpereV)^(-1);
% %ww = -ones(size(w))*max(w)-w;
% interpselfE = griddedInterpolant(w,selfE,'linear');
% dE = [-fliplr(w),w];
% % spacing = w(2)-w(1);
% for wi = w
%     i = i+1;
% %     integrand1 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
% %         +1).^(-1))/(wi+1i*impScat);
%     
%     k = 0; %initiate index corresponding to selfE(k) where k corresponds to epsilon
%     integrand2 = zeros(1,length(dE)); %initiate integrand for nonzero selfE epsilon>0
%     for ep = dE
%         k = k+1;
%         
%         if (ep+wi<0)
% %             if (ep<0)
%             integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%                 +1)^(-1))/(wi+interpselfE(-(wi+ep))-conj(interpselfE(-ep))+1i*impScat);
% %             else 
% %             integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
% %                 +1)^(-1))/(wi+interpselfE(-(wi+ep))+conj(interpselfE(ep))+1i*impScat);
% %             end
%         else
%             if (ep<0)
%             integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%                 +1)^(-1))/(wi-interpselfE(wi+ep)-conj(interpselfE(-ep))+1i*impScat);
%             else
%             integrand2(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
%                 +1)^(-1))/(wi-interpselfE(wi+ep)+conj(interpselfE(ep))+1i*impScat);
%             end
%         end
%     end
% %     integrand6 = zeros(1,length(w));
% %     k = 0;
% %     for epsi = ww
% %         k = k+1;
% %         integrand6(k) = ((exp(beta*(wi+epsi))+1).^(-1)-(exp(beta*epsi) ...
% %         +1).^(-1))/(wi+interpselfE(-(wi+epsi))+1i*impScat);
% %     end 
% % 
% %     integrand3 = zeros(1,length(w));
% %     k = 0;
% %     for eps = w
% %         k = k+1;
% %         
% %         if (eps+wi<=max(w)) 
% %             integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
% %         +1)^(-1))/(wi-interpselfE(wi+eps)+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e<=wmax
% %         
% %         else
% %             integrand3(k) = ((exp(beta*(wi+eps))+1)^(-1)-(exp(beta*eps) ...
% %         +1)^(-1))/(wi+conj(selfE(k))+1i*impScat); %defines integrand for selfE not equal to 0 w+e>wmax
% %         end
% %     end
% 
% %     integrand4 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
% %         +1).^(-1))/(wi-selfE(i)-imag(selfE(1))+1i*impScat);
% %     integrand5 = @(epsilon) ((exp(beta*(wi+epsilon))+1).^(-1)-(exp(beta*epsilon) ...
% %         +1).^(-1))/(wi+selfE(i)+imag(selfE(1))+1i*impScat);
% %     integrand3 = zeros(1,i-1); %initiate integrand for nonzero selfE epsilon<0
% %     
% %     k=0;
% %     for ep = -w(1:i-1)
% %         k=k+1;
% %         integrand3(k) = ((exp(beta*(wi+ep))+1)^(-1)-(exp(beta*ep) ...
% %         +1)^(-1))/(wi+selfE(-k+i)+1i*impScat);
% %     end
%     
%     cond(i) = kk*wp^2/(4i*pi*wi)*(trapz(dE,integrand2)); % + ... %Calculation of cond as a function of w,T for piecewise integral (selfE nonzero epsilon>0)
%         %integral(integrand4,0,min(w)) + ...
%         %integral(integrand5,-min(w),0)); %selfE zero    
%         %... %trapz(ww,integrand6)); % + ...
%         %trapz(w,integrand3)); % + ... %); % + ... %selfE nonzero epsilon<0
%         %integral(integrand1,-inf,-2*max(w)) + ... %selfE zero
%         %integral(integrand1,max(w),inf)); % + ... % + ...
%         %... %integral(integrand1,-min(w),min(w)) +...
% 
%     
% 
% end
% end
% 

