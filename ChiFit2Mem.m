function mem = ChiFit2Mem(x,rpaChi,Chi,w,kernel,T,dstart,dend)
%PitoDielectric Converts a bosonic spectral density into a dielectric function of w (eV)
%               for given values of w, and at temperature T (K) with plasma
%               frequency (wp) and impurity scattering rate (impScat)
%   Converts Pi to selfE, selfE to cond, cond to dielec
impScat = 0;
einf = 0;
wp = x(3);

% osc = [4264,22041800,4069;6789,8142014,3925;11650,5307060,3500;15409,45542000,8905;21300,223159025,13898;30756,320536896,6908;34946,214474009,6396;40421,756371984,7518];
% osc = osc/8065.6;
% osc(:,2) = osc(:,2)/8065.6;

% oscE = [0.48,1.34,0.9;1.41,0.44,0.67;2.51,1.54,1.78;3.36,0.23,0.24;4.87,3.92,2.21];
osc = [1.28320465902843,1.73993717599773,2.34492909044301,2.64845782511498;2.46325602296711,1.49809047299572,1.66488020974720,2.64845782511498;3.81921295116640,1.72246532204933,0.778691572537650,2.64845782511498;4.94333575701083,4.99990924685092,2.00918704525150,2.64845782511498];

epsL = LorentzOscill(w,osc(:,1),osc(:,2),osc(:,3),osc(1,4));

selfE = x(1)*KernelPitoSelfE(rpaChi,kernel,w)+x(2)*KernelPitoSelfE(Chi,kernel,w);
cond = SelfEtoCondNegative(selfE,w,T,impScat,wp);
dielec = CondtoDielec(cond,w,einf)+epsL;
cond = DielectoCond(w,dielec);
mem = CondtoMem(cond,w);
mem = mem(dstart:dend);
end