function [Pi] = GenerateBlocks(x)
%Block Generator generates a block(wprime) of width wii-wi and height fi
% wi and wii are in units of eV
% hbar = 1.05457173e-34; %planck's constant
% eV = 1.60217657e-19; %one eV in J
% wi = wi/hbar*eV;
% wii = wii/hbar*eV;

if length(x) == 11 %5 blocks
        Pi = @(wprime) (wprime>=x(1))*x(2) + (wprime>=x(3))*(-x(2)+x(4)) ...
                        + (wprime>=x(5))*(x(6)-x(4)) + (wprime>=x(7))*(x(8)-x(6)) ...
                        + (wprime>=x(9))*(x(10)-x(8)) - (wprime>=x(11))*x(10);    
elseif length(x) ==9 %4 blocks
        Pi = @(wprime) (wprime>=x(1))*x(2) + (wprime>=x(3))*(-x(2)+x(4)) ...
                        + (wprime>=x(5))*(x(6)-x(4)) + (wprime>=x(7))*(x(8)-x(6)) ...
                        - (wprime>=x(9))*x(8);
elseif length(x) == 7 %3 blocks
        Pi = @(wprime) (wprime>=x(1))*x(2) + (wprime>=x(3))*(-x(2)+x(4)) ...
                        + (wprime>=x(5))*(x(6)-x(4)) - (wprime>=x(7))*x(6);
elseif length(x) == 5% 2 blocks
                Pi = @(wprime) (wprime>=x(1))*x(2) + (wprime>=x(3))*(-x(2)+x(4)) ...
                        + (wprime>=x(5))*(-x(4));
elseif length(x) == 3% 1 block
                Pi = @(wprime) (wprime>=x(1))*x(2) - (wprime>=x(3))*x(2);
else
    error('myApp:argChk', 'Wrong number of input arguments')    
end

