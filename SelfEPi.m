function [selfE] = SelfEPi(kernel,Pi)
%SelfEnergyFunction Calculates the self energy from the kernel relating
%the bosonic spectrum to optical data. kernel is the Kernel Function
%relating bosonic excitations to ellipsometry and conductivity data. Pi is
%the bosonic glue function.

%X(w)
dw = 1e-3;
x1 = dw:dw:0.5; %w' (0,0.5)

%Initiate
selfE = zeros(4001,1); %selfE(w) w = (-1,1)

%Terms
Piwp = Pi(x1);

for i = 1:4001
    selfE(i) = dw*trapz(kernel(i,:).*Piwp);
end

end

