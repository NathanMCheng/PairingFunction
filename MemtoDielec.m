function [d] = MemtoDielec(mem,w,wp,einf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

d = einf-wp^2./(w.*(mem+w));

end

