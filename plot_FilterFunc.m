% script name: "plot_FilterFunc"
% Plotting the filterfunction
nameit = 'filterFunction';
c1 = 0.1;
r1 = .1;
R1 = .1;

c2 = 0.4;
r2 = .1;
R2 = .2;

fx1 = @(x) FilterFunc(x,R1,r1,c1);
fx2 = @(x) FilterFunc(x,R2,r2,c2);
plot_dom = -1:0.001:1;

figure;
h1 = plot(plot_dom,fx1(plot_dom),'k','LineWidth',4);
hold on;
h2 = plot(plot_dom,fx2(plot_dom),'-.r','LineWidth',5);
func1Name = num2str([c1, r1, R1] );  %'(0,0.3,0.3)';
func2Name = num2str([c2, r2, R2] );
h = legend(func1Name,func2Name,'Location','NorthWest','Interpreter','latex')
title(h, {'The parameters:' ; '(c)enter    (r)ise     (R)adius'})
%h = legend([h1, h2], func1Name,func2Name,'Location','northoutside','Orientation','horizontal');
set(gcf, 'Position', get(0,'ScreenSize'));
set(gca,'FontSize',27);
h.FontSize = 32;
grid on

saveas(gcf,nameit,'fig');        
saveas(gcf,nameit,'pdf');
print('-depsc2',nameit);