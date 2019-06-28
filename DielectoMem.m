function mem = DielectoMem(e1,e2,w,wp,einf)
%ConductivityToMemoryFunction 
% wp = 1; do not calculate with it

mem1 = wp^2./w.*(einf-e1)./((einf-e1).^2+e2.^2)-w;
mem2 = wp^2./w.*e2./((einf-e1).^2+e2.^2);

% ed = e1+1i*e2-einf+1;
% 
% if iscolumn(w)
%     w = w';
% end
% if iscolumn(ed)
%     ed = transpose(ed);
% end
% we = w.*(-ed);
% mem = (wp^2./we-w);
mem = mem1+1i*mem2;
% mem1 = real(mem);
% mem2 = imag(mem);
end

