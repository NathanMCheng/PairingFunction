function [Kp,Km] = GapKernel0(Pi)
%GapKernelFunction Kernel Function of the Gap Eq. at 0K
%   See PhysRev.148.263
ipi = 1i*pi();
dw = 1e-3;
v = dw:dw:0.4;
alphaF = Pi(v);

x1 = dw:dw:8;
x2 = -4+dw:dw:4;

K1 = zeros(8000,1);
K2 = zeros(8000,1);

Kp = zeros(4000,4000);
Km = zeros(4000,4000);

%no singularities
tic
for j = 4000:8000
    K1(j) = dw*trapz(alphaF./(x1(j)+v));
    K2(j) = dw*trapz(alphaF./(x2(j)+v));
end
toc

tic
for j = 1:3599
    K1(j) = dw*trapz(alphaF./(x1(j)+v));
    K2(j) = dw*trapz(alphaF./(x2(j)+v));
end
toc

%singularities
tic
%x2(3600) = -0.4
K1(3600) = dw*trapz(alphaF./(x1(3600)+v));
e = dw;
K0 = dw*trapz(alphaF(1:399)./(x2(3600)+v(1:399)))+ipi*alphaF(400);
while abs(K2(3600)-K0) > 1e-4
    K2(3600) = K0;
    e = e/10;
    eps1 = (v(400)-9*e):e:(v(400)-e);
    eps2 = (v(400)+e):e:(v(400)+9*e);
    K0 = K0+e*trapz(Pi(eps1)./(x2(3600)+eps1))+ ...
        e*trapz(Pi(eps2)./(x2(3600)+eps2));
end
K2(3600) = K0;

%x2(3999) = -dw
K1(3999) = dw*trapz(alphaF./(x1(3999)+v));
e = dw;
K0 = dw*trapz(alphaF(2:400)./(x2(3999)+v(2:400)))+ipi*alphaF(1);
while abs(K2(3999)-K0) > 1e-4
    K2(3999) = K0;
    e = e/10;
    eps1 = (v(1)-9*e):e:(v(1)-e);
    eps2 = (v(1)+e):e:(v(1)+9*e);
    K0 = K0+e*trapz(Pi(eps1)./(x2(3999)+eps1))+ ...
        e*trapz(Pi(eps2)./(x2(3999)+eps2));
end
K2(3999) = K0;

for j = 3601:3998 %x
    e = dw;
    i_0 = 4000-j;
    i_1 = i_0-1;
    i_2 = i_0+1;
    K1(j) = dw*trapz(alphaF./(x1(j)+v));
    denom2(1:i_1) = alphaF(1:i_1)./(x2(j)+v(1:i_1));
    denom2(i_2:400) = alphaF(i_2:400)./(x2(j)+v(i_2:400));
    
    K0 = dw*trapz(denom2)+ipi*alphaF(i_0);
    while abs(K2(j)-K0) > 1e-4
        K2(j) = K0;
        e = e/10;
        eps1 = (v(i_0)-9*e):e:(v(i_0)-e);
        eps2 = (v(i_0)+e):e:(v(i_0)+9*e);
        K0 = K0+e*trapz(Pi(eps1)./(x2(j)+eps1))+ ...
            e*trapz(Pi(eps2)./(x2(j)+eps2));
    end
    K2(j) = K0;
end
toc

%calculate Km, Kp
tic
for i = 1:4000 % w
    for j = 1:4000 %wp
        Kp(i,j) = K1(i+j)+K2(4000-i+j);
        Km(i,j) = K1(i+j)-K2(4000-i+j);
    end
end
toc

end

