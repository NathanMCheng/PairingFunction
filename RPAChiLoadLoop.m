function [rpaChi] = RPAChiLoadLoop(prefix,loops)
%RPAChiLoadLoop 

rpaChi = zeros(1001,1001,loops);

for i = 1:loops
    filename = strcat(prefix,int2str(i-1));
    rpaChi(:,:,i) = dlmread(filename,'\t',[0 0 1000 1000]);

end

