function [finalChi] = MomentumIntegrate(momentumDependent,x,y)
%MomentumIntegration 

%l = size(momentumDependent(:,1,1));
%w = size(momentumDependent(1,:,1);
h = length(momentumDependent(1,1,:));

%intermediateChi = zeros(1,l);
finalChi = zeros(1,h);

for i = 1:h
    intermediateChi = trapz(y,momentumDependent(:,:,i));
    finalChi(1,i) = trapz(x, intermediateChi);
end


end

