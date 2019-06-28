function [dtild,wtild,flag,iter] = elireal(Pi,T,iwtild,idtild,wtild0,dtild0)
%Real Axis Eliashberg
%   Analytic continuation of the imaginary axis Eliashberg equations.
flag = 0;
iwtild = iwtild(1002:2003);
idtild = idtild(1002:2003,:);
wtild = zeros(801,1);
dtild = zeros(801,314);

g=1;
KpereV = 11604.505;
kbt = T/KpereV;
pt = pi()*kbt;
delta = 1e-3;

m = 0:1001; %fermi frequency index +/-
ff = 1i*pt*(2*m+1); %fermi frequency

n = -1001:1001;
wm = 1i*pt*(2*n+1);
wn = wm;
clear m

dthet = 1e-2; %theta
thet = 0:dthet:pi()-dthet;

dz = 1e-3; %z
z = -0.4:dz:0.4;
fzv = zeros(801,801);
omg = 0:dz:0.4;

wtild0i = griddedInterpolant(z,wtild0,'linear','nearest');
[dtildgrid1,dtildgrid2] = ndgrid(z,thet);
dtild0i = griddedInterpolant(dtildgrid1,dtildgrid2,dtild0,'linear','nearest');

queryp = -0.8:1e-3:0.8;
[dtildgridq1,dtildgridq2] = ndgrid(queryp,thet);

dtild1 = dtild0i(dtildgridq1,dtildgridq2);
wtild1 = wtild0i(queryp);
clear wtild0i
clear dtild0i

% 'note M should be negative/positive here';
idtildavg = zeros(1002,1); %imaginary dtild avg
iwtildavg = zeros(1002,1); %imaginary wtild avg
idenom = zeros(1002,314); %imaginary wtild/dtild avg denom
lmbdmn = zeros(1002,2003); %imaginary lambda(m-n)

ipiomg = Pi(omg); %imaginiary I^2chi

% 'note M is only positive here
dtildavg = zeros(801,801); % dtild avg
wtildavg = zeros(801,801); % wtild avg
wtilddenom = zeros(801,801,314); % wtild avg denom
dtilddenom = zeros(801,801,314); % dtild avg denom

lmbdp = zeros(801,1002);
lmbdm = zeros(801,1002);

ichifb = zeros(801,801);

wtildsum = zeros(801,1);
dtildsum = zeros(801,1);

cos2thet = cos(2*thet);
piomg = Pi(z);
piomg(1:400) = -flipud(piomg(402:801));
piomg(401) = 0;


% f(z-v)
for i = 1:801 % z index
    fzv(:,i) = 1./(exp((z(i)-z)/kbt)+1);
end

% n(z)
nz = 1./(exp(z/kbt)-1);

% I^2*Chi[n(z)+f(z-v)] calculation
for i = 1:801 % v index
    ichifb(i,:) = piomg.*(nz+fzv(i,:));
end

ichifb(:,401) = 0;

% imaginary functions denominator <>'
for i = 1:1002 % wm index
    idenom(i,:) = sqrt(iwtild(i)^2+idtild(i,:).^2); %plus in all 4 papers
end

% imaginary wtild & dtild angle averages <>'
for i = 1:1002 % wm index
    idtildavg(i) = mean(idtild(i,:).*cos(2*thet)./idenom(i,:));
    iwtildavg(i) = mean(iwtild(i)./idenom(i,:));
end

% lambda m-n calculation m = -1001:1001 n = 0:1001
for i = 1:1002 % wn index
    for j = 1:2003 % wm index
        lmbdmnint = ipiomg.*omg./(omg.^2+(

clear idenom
clear thet
clear idtild
clear iwtild

% lmbd(v-iwmn)
for i = 1:801 % v index
    for j = 1:1002 % wm index
        lmbdint = piomg./(z(i)-ff(j)-z+1i*delta);
        lmbdm(i,j) = dz*trapz(lmbdint);
    end
end

% lmbd(v+iwmn)
for i = 1:801 % v index
    for j = 1:1002 % wm index
        lmbdint = piomg./(z(i)+ff(j)-z+1i*delta);
        lmbdp(i,j) = dz*trapz(lmbdint);
    end
end

clear lmbdint
clear piomg
clear ff

% wtild sum calculation
for i = 1:801 % v index
    wtildsum(i) = sum((lmbdm(i,:)-lmbdp(i,:)).*transpose(iwtildavg));
end

clear iwtildavg

% dtild sum calculation
for i = 1:801 % v index
    dtildsum(i) = sum((lmbdm(i,:)+lmbdp(i,:)).*transpose(idtildavg));
end

clear idtildavg
clear lmbdm
clear lmbdp

% wtild(v+id) denom calculation
for i = 1:801 % v index
    for j = 1:801 % z index
            wtilddenom(i,j,:) = ...
                sqrt(wtild1(801+i-j)^2-dtild1(801+i-j,:).^2);
    end
end

% wtildavg calculation
for i = 1:801 % v index
    for j = 1:801 % z index
        wtildavg(i,j) = mean(wtild1(801+i-j)./wtilddenom(i,j,:));
    end
end

% wtild calculation
for i = 1:801 % v index
    wtildint = 1e-3*trapz(ichifb(i,:).*wtildavg(i,:));
    wtild(i) = z(i)+1i*pt*wtildsum(i)+1i*pi()*wtildint;
end

wtildi = griddedInterpolant(z,wtild,'linear','nearest');
wtild1 = wtildi(queryp);
clear wtildi

% dtild(v+id) denom calculation
for i = 1:801 % v index
    for j = 1:801 % z index
        dtilddenom(i,j,:) = ...
            sqrt(wtild1(801+i-j)^2-dtild1(801+i-j,:).^2);
    end
end
clear wtild1

% dtildavg calculation
for i = 1:801 % v index
    for j = 1:801 % z index
        dtildavg(i,j) = ...
            mean(transpose(dtild1(801+i-j,:).*cos2thet)./squeeze(dtilddenom(i,j,:)));
    end
end
clear dtild1

% dtild calculation
for i = 1:801 %v index
    dtildint = cos2thet.*1e-3*trapz(ichifb(i,:).*dtildavg(i,:));
    dtild(i,:) = 1i*pt*g*dtildsum(i).*cos2thet+1i*pi()*g*dtildint;
end

iter = 1;
while norm(abs(dtild-dtild0)) > 1e-6
    
    if iter>=100 % max iter reached
        flag = 1;
        break
    elseif norm(dtild) <1e-5 % dtild converged to 0
        flag = 2;
        break
    elseif norm(wtild) <1e-5 % wtild convered to 0
        flag = 3;
        break
    end
    
    % update dtild0 & wtild0, iter
    wtild0 = wtild;
    dtild0 = dtild;
    iter = iter+1;
    
    wtild0i = griddedInterpolant(z,wtild0,'linear','nearest');
    dtild0i = griddedInterpolant(dtildgrid1,dtildgrid2,dtild0,'linear','nearest');
    
    dtild1 = dtild0i(dtildgridq1,dtildgridq2);
    wtild1 = wtild0i(queryp);
    clear wtild0i
    clear dtild0i
        
    % wtild(v+id) denom calculation
    for i = 1:801 % v index
        for j = 1:801 % z index
            wtilddenom(i,j,:) = ...
                sqrt(wtild1(801+i-j)^2-dtild1(801+i-j,:).^2);
        end
    end
    
    % wtildavg calculation
    for i = 1:801 % v index
        for j = 1:801 % z index
            wtildavg(i,j) = mean(wtild1(801+i-j)./wtilddenom(i,j,:));
        end
    end
    
    % wtild calculation
    for i = 1:801 % v index
        wtildint = 1e-3*trapz(ichifb(i,:).*wtildavg(i,:));
        wtild(i) = z(i)+1i*pt*wtildsum(i)+1i*pi()*wtildint;
    end
    
    wtildi = griddedInterpolant(z,wtild,'linear','nearest');
    wtild1 = wtildi(queryp);
    clear wtildi
    
    % dtild(v+id) denom calculation
    for i = 1:801 % v index
        for j = 1:801 % z index
                dtilddenom(i,j,:) = ...
                    sqrt(wtild1(801+i-j)^2-dtild1(801+i-j,:).^2);
        end
    end
    clear wtild1
    
    % dtildavg calculation
    for i = 1:801 % v index
        for j = 1:801 % z index
            dtildavg(i,j) = ...
                mean(transpose(dtild1(801+i-j,:).*cos2thet)./squeeze(dtilddenom(i,j,:)));
        end
    end
    clear dtild1
    
    % dtild calculation
    for i = 1:801 %v index
        dtildint = cos2thet.*1e-3*trapz(ichifb(i,:).*dtildavg(i,:));
        dtild(i,:) = 1i*pt*g*dtildsum(i).*cos2thet+1i*pi()*g*dtildint;
    end
    
end

