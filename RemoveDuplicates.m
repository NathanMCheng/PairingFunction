function [xprime,yprime] = RemoveDuplicates(x,y)
%RemoveDataDuplicates
lngth = length(x);
xprime = x;
yprime = y;
k = 1;
while k < lngth
    if round(1000*xprime(k))/1000 == round(1000*xprime(k+1))/1000
        lngth = lngth-1;
        xprime(k) = [];
        yprime(k) = [];
    end
    k = k+1;
end

        
end

