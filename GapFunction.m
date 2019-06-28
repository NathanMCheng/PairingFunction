function [delta,Z,delta0,flag,iter,step] = GapFunction(delta0,Kp,Km,step)
%GapFunction Calculates Delta and Z from kernels
%   see PhysRev.148.263
flag = 0;

dw = 1e-3;
w = dw:dw:4;
w = w';
wp = dw:dw:4;
wp = wp';

delta = zeros(4000,1);
Z = zeros(4000,1);
tic
% Z and Delta Calculation
denom = (wp.^2-delta0.^2).^(1/2);
ReZ = real(wp./denom);
ReD = real(delta0./denom);
for i = 1:4000
    Z(i) = 1-dw*trapz(ReZ.*transpose(Km(i,:)))/w(i);
    delta(i) = 1/Z(i)*dw*trapz(ReD.*transpose(Kp(i,:)));
end
toc
iter = 1;
norm1 = norm(abs(delta-delta0));
norm0 = norm1;
% step = 0.1;
while norm1 > 1e-5
    if iter == 100
        flag = 1;
        break;
    end
    if norm(abs(Z)) < 1e-3
        flag = 2;
        break;
    end
    if norm(abs(delta)) < 1e-5
        flag = 3;
        break;
    end
        learn = norm1/norm0;
        if learn <= 1
            learn2 = 0;
            if step < 5e-5
                step = 1e-2;
            else
                step = step/1.25;
            end
            delta0 = delta0*(1-step)+delta*step;
        else
            learn2 = learn2+1;
            if learn2 >= 2
                step = 0.05;
                delta0 = delta0*(1-step)+delta*step;
    %             delta0 = delta;
            else
                delta0 = delta0*(1-step)+delta*step;
            end
        end
    delta0 = delta0*(1-step)+delta*step;
    iter = iter+1;
    
    % Z and Delta Calculation
    denom = (wp.^2-delta0.^2).^(1/2);
    ReZ = real(wp./denom);
    ReD = real(delta0./denom);
    for i = 1:4000
        Z(i) = 1-dw*trapz(ReZ.*transpose(Km(i,:)))/w(i);
        delta(i) = 1/Z(i)*dw*trapz(ReD.*transpose(Kp(i,:)));
    end

    norm0 = norm1;
    norm1 = norm(abs(delta-delta0));
    max(abs(delta-delta0));
    d0 = real(delta(1))
    x = (w-d0)/0.035;
    plot(x,real(delta)/d0,'b',x,(imag(delta))/d0,'r')
    axis([0,5,-inf,inf])
    drawnow
end


            




end

