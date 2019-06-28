function [sigma,w,theta] = ReftoCond(R,wprime)
% Converts measured reflectivity to conductivity data
%   Given reflectance (R) and corresponding freqeuency (w), approximates the
%   phase (theta) of the reflectivity via a finite Kramers-Kronig integral and then
%   solves for conductivity (sigma) given the relation between conductivity and
%   reflectivity


theta=zeros(size(wprime)); %initialize phase
sigma=zeros(size(wprime)); %initialize conductivity

i=0; %initialize index of sigma and theta. determines correspondence between frequency and conductivity/phase
for w = wprime %loop for all values of frequency
    i=i+1;
    f = (log(R)-log(R(i)))./(wprime.^2-w^2);
    
    for k = 1:length(wprime) %identify singularities
        if isnan(f(k))
            f(k)=0;
        end
    end
        
    theta(i) = -w/pi()*trapz(wprime,f); %integrate and solve for reflectivity phase
    r = R(i)^(1/2)*exp(1i*theta(i)); %determine the reflectivity coefficient from the reflectance and phase
    sigma(i) = 1i*w/4/pi()*(1-((1+r)/(r-1))^2); %solve for conductivity
        
end 
w=wprime;
end

