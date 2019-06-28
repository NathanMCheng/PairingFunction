function [mkq] = QCoupling(t)
%mkqCoupling Calculates Q variant coupling according to Keimer's paper
%

%Variables with given values
J = .13; %eV
%wqAF = 0.06; %eV
%GammaQ = 0.1; %eV
%

Qy = -pi():0.1:pi(); %Qy momentum
Qx = -pi():0.1:pi(); %Qx momentum

kY = -pi():0.1:pi(); %kY momentum
kX = -pi():0.1:pi(); %kX momentum

a = 0; %Qy index
b = 0; %Qx index
c = 0; %kX index
%d = 0; %kY index

%gammaQ = zeros(length(Qx),length(Qy));
%ek = zeros(length(kX),length(kY));
%ekq = zeros(length(kX),length(kY),length(Qx),length(Qy));
mkq = zeros(length(kX),length(kY),length(Qx),length(Qy));

for qx = Qx
    b  = b+1; %Qx index
    
    for qy = Qy
        a = a+1;
        gammaQ = (cos(qx)+cos(qy))/2; 
        
        for ki = kX
            c = c+1;
            ek = -2*t*(cos(ki)+cos(kY))-4*(-t/3)*cos(ki)*cos(kY);
            kqx = ki-qx; %kx-qx
            kqY = kY-qy; %ky-qy
            ekq = -2*t*(cos(kqx)+cos(kqY))-4*(-t/3)*cos(kqx)*cos(kqY);
            mkq(a,b,c,:) = 2*J*gammaQ*ones(size(ek)) + 1/2*(ekq+ek);
        end
    end
end






end

