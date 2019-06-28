function [deltaR] = FittodeltaR(w,R0,R1,Noneqkernel,Pibe,Piscp,Pilat,Te,Tbe,Tscp,Tlat,wR)
%FitToReflectivityChange Fit to the change in reflectivity using the 4TM
%model at a given frequency
%


impScat = 0;
wp = 1;
selfE = NoneqKernelPitoSelfE(Pibe,Noneqkernel,w,Te,Tbe) + ...
    NoneqKernelPitoSelfE(Piscp,Noneqkernel,w,Te,Tscp) + ...
    NoneqKernelPitoSelfE(Pilat,Noneqkernel,w,Te,Tlat);
cond = NoneqSelfEtoCondNegative(selfE,w,Te,impScat,wp,wR);
Reflec = transpose(CondtoRef(cond,wR));
deltaR = (Reflec-R0)/R1;

end

