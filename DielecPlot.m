function [] = DielecPlot(f,T,sample,w,e1,e2,fitw,fit)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

figure;
fig = plot(w,e1,'*b');
hold on
fig2 = plot(w,e2,'*b');
p1 = plot(fitw,real(fit),'r','LineWidth',1.5);
p2 = plot(fitw,imag(fit),'r','LineWidth',1.5);
axis([0,1,-100,100])
xlabel('Energy (eV)')
ylabel('Dielectric Function')
legend([fig,p1],sample,[f ' Fit'])

saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f T '.eps'],'epsc');
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f T '.fig']);
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f T '.jpg']);

end

