function [Chi] = RPAChiNormal(frequency,qx,qy)
%RENORMALIZED-BCS-LINHARDFUNCTION Renormalized random phase approximation
%in terms of the BCS-Lindhard response function
%   Chi '' taken from Thesis: Andreas Schnyder
%
%   K is a variable existing in the Brillouin zone: K(1) = Kx, K(2) = Ky.
%   Q is a scattering variables: Q(1) = Qx, Q(2) = Qy
%

%table variables
% delta0 = 0.035; %eV (35.0 meV) superconducting gap
chemPotent = 0.0879; %eV (87.9 meV)
t1 = -0.5547; %eV (-554.7 meV)
t2 = 0.1327; %eV (132.7 meV)
t3 = 0.0132; %eV (13.2 meV)
t4 = -0.1849; %eV (184.9 meV)
t5 = 0.0265; %eV (26.5 meV)
U = 0.165; %eV (165.0 meV) U-limit of the Hubbard model

%additional variables needed
KpereV = 11604.505;
beta = 1/(300/KpereV); %1/eV
Gamma = 0.0024; %eV (1 meV)

%k and q variables
kx = -pi():2*pi()/1000:pi();
ky = -pi():2*pi()/1000:pi();
% qx = 0.5*pi():pi()/50:1.5*pi();
% qy = 0.5*pi():pi()/50:1.5*pi();
N = length(kx)*length(ky); % number of sites on a square lattice wiht periodic boundary conditions

Chi = zeros(length(frequency),length(qx),length(qy)); % Chi'' initation
l = 0; %w index
%m = 0; %qx index
%n = 0; %qy index
for w = frequency
    l = l+1;
    m = 0; % qx index
     for qxi = qx
        m = m+1;
        n = 0; %qy index
        for qyi = qy
            n = n+1;
            sum = 0;
            for kxi = kx
                for kyi = ky
                    
                    cosKx = cos(kxi);
                    cosKy = cos(kyi);
                    cos2Kx = cos(2*kxi);
                    cos2Ky = cos(2*kyi);
                    
                    cosKQx = cos(qxi+kxi);
                    cosKQy = cos(qyi+kyi);
                    cos2KQx = cos(2*(qxi+kxi));
                    cos2KQy = cos(2*(qyi+kyi));
                    
                    epsilonK =  t1/2*(cosKx+cosKy)+t2*cosKx*cosKy +t3/2*(cos2Kx+cos2Ky) + ...
                        t4/2*(cos2Kx*cosKy+cosKx*cos2Ky)+t5*cos2Kx*cos2Ky+chemPotent;
                    %band parameters at various points in the BZ
                    
%                     deltaK = delta0/2*(cosKx-cosKy); %gap at various points in the BZ
                    
                    epsilonKQ =  t1/2*(cosKQx+cosKQy)+t2*cosKQx*cosKQy +t3/2*(cos2KQx+cos2KQy) + ...
                        t4/2*(cos2KQx*cosKQy+cosKQx*cos2KQy)+t5*cos2KQx*cos2KQy+chemPotent;
                    %band parameters at various points in the BZ
                    
%                     deltaKQ = delta0/2*(cosKQx-cosKQy); %gap at various points in the BZ
                    
%                     EK = sqrt(epsilonK^2+deltaK^2); %BCS dispersion of quasiparticles
                    
%                     EKQ = sqrt(epsilonKQ^2+deltaKQ^2); %BCS dispersion of qusiparticles
                    
%                     CplusQK = (1 + (epsilonKQ*epsilonK+deltaKQ*deltaK)/(EKQ*EK)); %coherence factor
                    
%                     CminusQK = (1 - (epsilonKQ*epsilonK+deltaKQ*deltaK)/(EKQ*EK)); %coherence factor
                    
                    fermiEKQ = 1/(exp(beta*epsilonKQ)+1);
                    
                    fermiEK = 1/(exp(beta*epsilonK)+1);
                    
%                     fermiMinusEK = 1/(exp(-beta*EK)+1);
                    
%                     fermiMinusEKQ = 1/(exp(-beta*EKQ)+1);
                    
                    
                    sum =  sum + (fermiEKQ-fermiEK)/(w+1i*Gamma-(epsilonKQ-epsilonK));
%                     sum = sum + 1/2*CplusQK*(fermiEKQ-fermiEK)/(w+1i*Gamma-(EKQ-EK)) ...
%                         + 1/4*CminusQK*(1-fermiEKQ-fermiEK)/(w+1i*Gamma+(EKQ+EK)) ... 
%                         + 1/4*CminusQK*(fermiEKQ+fermiEK-1)/(w+1i*Gamma-(EKQ+EK));
                end
            end
            
            chi0 = sum/N; %Chi0
%             Chi(l,m,n) = imag(chi0/(1-U*chi0));
            Chi(l,m,n) = imag(chi0)/((1-U*real(chi0))^2+U^2*imag(chi0)^2); %Chi''
        end
    end
end

end
