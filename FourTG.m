function [G] = FourTG(Pi,w,Tb,Te)
%FourTemperatureG Calculates G of the four temperature model
%   Detailed explanation goes here
gammae = 1e-4; % J/cm^3/K^2
kb = 8.61733e-5; % eV/K
hbar = 6.582119e-16; % eV*s
KpereV = 11604.505; % K/eV
integral = Pi(w).*w.^2.*((exp(w*KpereV/Tb)-1).^(-1)-(exp(w*KpereV/Te)-1).^(-1));
G = 6*gammae/(pi()*hbar*kb^2)*trapz(w,integral); % J/cm^3/s

end

