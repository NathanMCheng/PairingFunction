function [chi] = MomentumDependentChi(w,A)
%MomentumDependentChiCalculation Calculates the mommentum dependent 3-D Chi
%for various momentums and frequencies
% 

%Variables with given values
J = .13; %eV
wqAF = 0.06; %eV
GammaQ = 0.1; %eV
%

q = 0:0.02:pi(); %q momentum
k = 0:0.02:pi(); %k momentum

m = 0; %q index
n = 0; %k index

gammaQ = zeros(length(k),length(q));
for kn = k
    n  = n+1;
    gammaQ(n,:) = (cos(kn)+cos(q))/2; 
end
 
nq = A*(1-gammaQ);

wq = zeros(length(k),length(q));
n = 0; %reset k index
for kn = k
    n = n+1;
    wq(n,:) = 2*J*sqrt((1-gammaQ(n,:)).*(1+gammaQ(n,:)+wqAF^2/(8*J^2)));
end

chi = zeros(length(w),length(k),length(q));

for qm = q
    m = m+1;
    n = 0; %reset k index
    for kn = k
        n = n+1;
        chi(:,n,m) = nq(n,m)./(wq(n,m)^2-w.^2-1i*GammaQ*w);
    end
end
    







end

