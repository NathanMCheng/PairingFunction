function [] = ModelPlot(f,sample,w,fit)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

fig = figure;
plot(w,fit(w),'r','LineWidth',1.5);
axis([0,0.4,0,3])
xlabel('Energy (eV)','FontSize',20)
ylabel('\Pi(\omega)','FontSize',20)
ax1 = gca;
set(ax1,'OuterPosition',[0,0.025,1,1])
leg1 = legend(sample);
set(leg1,'FontSize',20)
set(ax1,'FontSize',20)
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f 'Pi.eps'],'epsc');
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f 'Pi.jpg']);
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f 'Pi.fig']);

end

