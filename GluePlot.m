function [] = GluePlot(f,T,sample,w,fit)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

fig = figure;
plot(w,fit(w),'r','LineWidth',1.5);
ax = gca;
axis([0,0.5,0,2])
xlabel('Energy (eV)','FontSize',16)
ylabel('\Pi(\omega)','FontSize',16)
leg = legend(sample);
set(leg,'FontSize',16)
set(ax,'FontSize',16)

saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f T 'Pi.jpg']);
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f T 'Pi.eps'],'epsc');
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f T 'Pi.fig']);


end

