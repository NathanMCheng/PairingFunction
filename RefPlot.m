function [] = RefPlot(f,sample,w,R,fit)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

fig = figure;
p1 = plot(w(1:359),R(1:359),'*b');
hold on
p2 = plot(w(1:359),fit,'r','LineWidth',1.5)
ax = gca;
axis([0,1,0,1])
xlabel('Energy (eV)','FontSize',20)
ylabel('Reflectivity','FontSize',20)
set(ax,'OuterPosition',[0,0.025,1,1],'FontSize',20)
leg = legend(sample,[f ' Fit']);
set(leg,'FontSize',20)

saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f '.eps'],'epsc');
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f '.fig']);
saveas(fig,['C:\Users\Nathan\Documents\MATLAB\Thesis\EPS Figs\' f '.jpg']);
end

