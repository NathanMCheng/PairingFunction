function [y] = mypsi(xx)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
    cof = [76.18009173D0,-86.50532033D0,24.01409822D0, -1.231739516D0,.120858003D-2,-.536382D-5];
    stp = 2.50662827465D0;
    
    x = xx-1;
    tmp = x+5.5;
    tmp = (x+1/2).*log(tmp)-tmp;
    ser = 1;
    for j = 1:6
        x = x+1;
        ser = ser+cof(j)./x;
    end
    y = tmp + log(stp.*ser);



end

