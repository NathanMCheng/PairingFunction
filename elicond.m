function [dtild,wtild] = elicond(dtild0,wtild0,Pi,T)
%GapFunction Calculates the gap funciton using coupled Eliashberg Eq.
%   Detailed explanation goes here

g=1;
KpereV = 11604.505;
T = T/KpereV;
pT = pi()*T;

m = -1001:1001;
ff = 1i*pT*(2*m+1);
% bf = 2i*pT*n;
% beta = (T)^(-1);
wf = imag(ff);

dw = 1e-3;
dthet = 1e-2;
thet = 0:dthet:pi()-dthet;
lmbOmg = dw:dw:0.4;

wtild0avg = zeros(2003,1);
dtild0avg = zeros(2003,1);
denom = zeros(2003,314);
denom2 = zeros(2003,314);

lbmn = zeros(2003,2003);
lmbin = lmbOmg.*Pi(lmbOmg);
lmbOmg2 = lmbOmg.^2;
% D = zeros(1001,1001);
% Omg = zeros(1001,1001);
% impt = zeros(1001,1001);
wtild = zeros(2003,1);
dtild = zeros(2003,314);

% ff convolution
fmn = zeros(2003,2003);
for n = 1:2003
    fmn(n,:) = (wf-wf(n)).^2;
end

% wtild0avg
for i = 1:2003
    denom(i,:) = sqrt(wtild0(i)^2+dtild0(i,:).^2);
    wtild0avg(i) = mean(wtild0(i)./denom(i,:));
end

% lbmn calculation
for k = 1:2003 %n
    for j = 1:2003 %m
        lmbi = lmbin./(lmbOmg2+fmn(j,k));
        lbmn(j,k) = 2*trapz(lmbOmg,lmbi);
%         D(j,k) = mean(dtild0(k,:)./denom(j,:));
%         Omg(j,k) = mean(wtild0(k)./denom(j,:));
%         impt(j,k) = pi()*Gmma*Omg(1,k)./(c^2+Omg(1,k).^2+D(1,k).^2);
    end
end

% wtild calculation
for k = 1:2003 %n
    wtild(k) = wf(k)+pT*sum(lbmn(:,k).*wtild0avg);
end

% dtild denom
for i = 1:2003
    denom2(i,:) = sqrt(wtild(i).^2+dtild0(i,:).^2);
    dtild0avg(i) = mean(cos(2*thet).*dtild0(i,:)./denom2(i,:));
end

% dtild calc
for k = 1:2003 %n 
    dtild(k,:) = g*pT*cos(2*thet)*sum(lbmn(:,k).*dtild0avg);
end