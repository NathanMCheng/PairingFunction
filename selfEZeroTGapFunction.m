function [dlt,z,flag,iter] = selfEZeroTGapFunction(dtild0,wtild0,Pi,T)
%GapFunction Calculates the gap funciton using coupled Eliashberg Eq.
%   Detailed explanation goes here
flag = 0;

KpereV = 11604.505;
T = T/KpereV;
pT = pi()*T;

m = 0:1000;
ff = 1i*pT*(2*m+1);
bf = 2i*pT*n;
beta = (T)^(-1);
wn = real(ff);



dw = 1e-3;
dthet = 1e-2;
thet = 0:dthet:pi()-dthet;
lmbOmg = 0:dw:0.4;

fmn = zeros(1001,1001);
for n = 0:1000
    fmn(n,:) = (ff(m)-ff(n)).^2;
end

lbmn = zeros(1001,1001);
lmbin = lmbOmg.Pi(lmbOmg);
lmbOmg2 = lmbOmg.^2;
D = zeros(1001,1001);
Omg = zeros(1001,1001);
impt = zeros(1001,1001);
wtild = zeros(1001,1/dthet);
denom = sqrt(wtild0.^2+dtild0.^2);
for k = 0:1000 %n
    for j = 0:1000 %m
        lmbi = lmbin./(lmbOmg2+fmn(j,k));
        lbmn(j,k) = 2*trapz(lmbOmg,lmbi);
        D(j,k) = mean(dtild0(k,:)./denom(j,:));
        Omg(j,k) = mean(wtild0(k)./denom(j,:));
        impt(j,k) = pi()*Gmma*Omg(1,k)./(c^2+Omg(1,k).^2+D(1,k).^2);
    end
    wtild(k,:) = wn(k)+pT*cos(2*thet)*sum(lbmn(:,k).*Omg(:,k));
end



% 
% d = 1e-5;
% mu = 0.0879;
% 
% wl = dw:dw:1;
% % fd = 1./(exp(beta*wl)+1);
% % fdn = 1./(exp(-beta*wl)+1);
% dlt0 = dlt0(wl);
% selfE = selfE(wl);
% z = (1+selfE)./wl;
% %indices
% l = 0; %w
% % n = 0; %v
% 
% v = 0:dw:0.4;
% % bd = 1./(exp(beta*v)-1);
% 
% i1 = zeros(1,1/dw);
% i2 = zeros(1,1/dw);
% % i3 = zeros(1,1/dw);
% dlt = zeros(1,1/dw); 
% % z = zeros(1,1/dw); 
% 
% for w = dw:dw:1
%     l = l+1;
%     m = 0; %wp
%     for wp = dw:dw:1
%         m = m+1;
%         i1(m) = real(dlt0(m)/sqrt(wp^2-dlt0(m)^2)); %int1
%         ii2 = pi(w)*(1./(w+wp+v+1i*d)-1./(w-wp-v+1i*d)- ...
%             1./(w+v-wp+1i*d)-1./(w-v+wp-1i*d));
%         i2(m) = trapz(v,ii2);
%         
% %         ii3 = pi(w)*(1./(w+wp+v+1i*d)+1./(w-wp-v+1i*d)+ ...
% %             1./(w+v-wp+1i*d)+1./(w-v+wp-1i*d));
% %         i3(m) = trapz(v,ii3);
% %         i4 = real(wp/sqrt(wp^2-dlt0(m)^2));
%     end
% %     z(l) = 1-(dw*trapz(i4.*i3))/w;
%     dlt(l) = (dw*trapz(i1.*i2)-mu*dw*trapz(i1))/z(l);
%  
% end
% iter = 1;
% while (norm(abs(dlt0-dlt))) > 1e-5 % && (norm((abs(z-z0)))) > 0.01
%     if iter == 10000
%         flag = 1;
%         break;
%     end
%     dlt0 = dlt;
% %     z0 = z;
%     iter = iter+1;
%     l = 0;
%     for w = dw:dw:1
%         l = l+1;
%         m = 0; %wp
%         for wp = dw:dw:1
%             m = m+1;
%             i1(m) = real(dlt0(m)/sqrt(wp^2-dlt0(m)^2)); %int1
%             ii2 = pi(w)*(1./(w+wp+v+1i*d)-1./(w-wp-v+1i*d)- ...
%                 1./(w+v-wp+1i*d)-1./(w-v+wp-1i*d));
%             i2(m) = trapz(v,ii2);
%             
% %             ii3 = pi(w)*(1./(w+wp+v+1i*d)+1./(w-wp-v+1i*d)+ ...
% %                 1./(w+v-wp+1i*d)+1./(w-v+wp-1i*d));
% %             i3(m) = trapz(v,ii3);
% %             i4 = real(wp/sqrt(wp^2-dlt0(m)^2));
%         end
% %         z(l) = 1-(dw*trapz(i4.*i3))/w;        
%         dlt(l) = (dw*trapz(i1.*i2)-mu*dw*trapz(i1))/z(l);
% 
%     end
% end
% 
end