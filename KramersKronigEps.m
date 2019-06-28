function [e2] = KramersKronigEps(e1,w,einf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
l = length(w);
e2 = zeros(1,l);
i = 0;
if length(einf)==1
    for wp = w
        i = i+1;
        if i==1 || i==2
            int2 = (e1(i+1:l)-einf)./(w(i+1:l)-wp);
            e2(i) = -1/pi()*(trapz(i+1:l,int2));
        elseif i ==l || i == (l-1)
            int1 = (e1(1:i-1)-einf)./(w(1:i-1)-wp);
            e2(i) = -1/pi()*(trapz(w(1:i-1),int1));
        else
            int1 = (e1(1:i-1)-einf)./(w(1:i-1)-wp);
            int2 = (e1(i+1:l)-einf)./(w(i+1:l)-wp);
            e2(i) = -1/pi()*(trapz(w(1:i-1),int1)+trapz(i+1:l,int2));
        end
    end
else
    for wp = w
        i = i+1;
        if i==1 || i==2
            int2 = (e1(i+1:l)-einf(i+1:l))./(w(i+1:l)-wp);
            e2(i) = -1/pi()*(trapz(i+1:l,int2));
        elseif i ==l || i == (l-1)
            int1 = (e1(1:i-1)-einf(1:i-1))./(w(1:i-1)-wp);
            e2(i) = -1/pi()*(trapz(w(1:i-1),int1));
        else
            int1 = (e1(1:i-1)-einf(1:i-1))./(w(1:i-1)-wp);
            int2 = (e1(i+1:l)-einf(i+1:l))./(w(i+1:l)-wp);
            e2(i) = -1/pi()*(trapz(w(1:i-1),int1)+trapz(i+1:l,int2));
        end
    end
end

