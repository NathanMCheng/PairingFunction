function [reflectivity,r] = DielectoRef(epsilon1,epsilon2)
%DielectrictoReflectivity Computes the Reflectivity(w) from the Dielctric(w)
dielec = epsilon1 + 1i*epsilon2;
r = (1-sqrt(dielec))./(1+sqrt(dielec));
reflectivity = r.*conj(r);

end

