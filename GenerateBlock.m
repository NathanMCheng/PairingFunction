function [block] = GenerateBlock(fi,wi,wii)
%Block Generator generates a block(wprime) of width wii-wi and height fi
% wi and wii are in units of eV
% hbar = 1.05457173e-34; %planck's constant
% eV = 1.60217657e-19; %one eV in J
% wi = wi/hbar*eV;
% wii = wii/hbar*eV;
block = @(wprime) (wprime>=wi)*fi - (wprime>wii)*fi; %define a single block of histogramic representation

end

