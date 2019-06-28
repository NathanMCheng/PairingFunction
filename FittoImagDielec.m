function [dielec,Pi] = FittoImagDielec(x,w,kernel)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
T = 300;
impScat = 0;
wp = 1;
if iscolumn(w)
    w = w';
end

if length(x) == 7
        Pi = @(wprime) (wprime>=x(1))*x(2) + (wprime>=x(3))*(-x(2)+x(4)) ...
                        + (wprime>=x(5))*(x(6)-x(4)) - (wprime>=x(7))*x(6); %define 3 blocks of histogramic representation
elseif length(x) == 5
                Pi = @(wprime) (wprime>=x(1))*x(2) + (wprime>=x(3))*(-x(2)+x(4)) ...
                        + (wprime>=x(5))*(-x(4));

elseif length(x) == 3
                Pi = @(wprime) (wprime>=x(1))*x(2) - (wprime>=x(3))*x(2);
else
    error('myApp:argChk', 'Wrong number of input arguments')    
end

selfE = KernelPitoSelfE(Pi,kernel,w);
cond = SelfEtoCond(selfE,w,T,impScat,wp);
dielec = imag(transpose(CondtoDielec(cond,w)));

end