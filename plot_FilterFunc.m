% script name: "plot_FilterFunc"
% Plotting the filterfunction
nameit = 'filterFunc';

fx1 = @(x) FilterFunc(x,.3,.2,0);
fx2 = @(x) FilterFunc(x,.2,.15,.5);
plot_dom = -1:0.001:1;

figure;
h1 = plot(plot_dom,fx1(plot_dom),'k','LineWidth',4);
hold on;
h2 = plot(plot_dom,func2(plot_dom),'-.r','LineWidth',5);
func1Name = 'c=0, R=.3, r=.2';
func2Name = 'c=.5, R=.2, r=.15';
h = legend([h1, h2], func1Name,func2Name,'Location','northoutside','Orientation','horizontal');
set(gcf, 'Position', get(0,'ScreenSize'));
set(gca,'FontSize',27);
grid on

saveas(gcf,nameit,'fig');        
saveas(gcf,nameit,'pdf');
print('-depsc2',nameit);