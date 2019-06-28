function [rpaChi] = LoadRPAChiText(filename)
%LoadFunctionforRPAChi

% rpaChi = zeros(blocksize,blocksize);
%sheet = 1;
% i = 0;
% for j = 1:blocks
%     i = i + 1;
%     if (i>blocks)
%         i = 1;
%         %sheet = sheet+1;
%
%     end
%sheet = strcat('Sheet',int2str(sheetindex));
%index = int2str((i-1)*401+1);
%     if (index == '0')
%         index = '1';
%     end
%index = ((i-1)*401+i);
%limit = int2str(i*401);
%limit = (i*401+(i-1));
%range = strcat('A',index,':OK',limit);

rpaChi(:,:) = dlmread(filename,'\t',[0 0 1000 1000]);

%rpaChi(:,:,j) = xlsread(filename,sheet,range);
end

