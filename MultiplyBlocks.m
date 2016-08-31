% script name: "MultiplyBlocks"
%
% We test effect of different Jordan blosck sizes
% The test functions: x./(x^2+1)

svait = 1;
name_it = 'multiBlocks';


f = @(x)  x./(x^2+1);

coef_num = 40;
cc    = chebcoefs(f,6000); 
coefs = cc(1:coef_num);

coefs(1:2:coef_num) = 0;


lambda = 0.5;

% Jordan blocks
J10 = lambda*eye(10)+.5*(triu(ones(10),1)-triu(ones(10),2));
J5 = lambda*eye(5)+.5*triu(ones(5),1)-triu(ones(5),2);
J2 = lambda*eye(2)+.5*triu(ones(2),1)-triu(ones(2),2);

% the matrices
A1 = J10;
A2 = blkdiag(J5,J5);
A3 = blkdiag(J2,J2,J2,J2,J2);
A4 = blkdiag(J10,J10);

cn_1 = cond(A1)
cn_2 = cond(A2)
cn_3 = cond(A3)
cn_4 = cond(A4)


f_A1 =  A1*(A1^2+eye(size(A1,1)))^(-1);
f_A2 =  A2*(A2^2+eye(size(A2,1)))^(-1);
f_A3 =  A3*(A3^2+eye(size(A3,1)))^(-1);
f_A4 =  A4*(A4^2+eye(size(A4,1)))^(-1);

number_of_plots = 30;
plot_dom = floor(linspace(5,coef_num,number_of_plots));

error_ar_1 = zeros(1,number_of_plots);
error_ar_2 = zeros(1,number_of_plots);
error_ar_3 = zeros(1,number_of_plots);
error_ar_4 = zeros(1,number_of_plots);

for l=1:number_of_plots
    Cheb = chebeval(coefs,A1,plot_dom(l));
    error_ar_1(l) = norm(Cheb-f_A1)/cn_1;
    Cheb = chebeval(coefs,A2,plot_dom(l));
    error_ar_2(l) = norm(Cheb-f_A2)/cn_2;
    Cheb = chebeval(coefs,A3,plot_dom(l));
    error_ar_3(l) = norm(Cheb-f_A3)/cn_3;
    Cheb = chebeval(coefs,A4,plot_dom(l));
    error_ar_4(l) = norm(Cheb-f_A4)/cn_4;
end;

figure();
h1 = semilogy(plot_dom,error_ar_1,'LineWidth',4);
hold on;
h2 = semilogy(plot_dom,error_ar_2,'--r','LineWidth',5);
h3 = semilogy(plot_dom,error_ar_3,':k','LineWidth',4);
%h4 = semilogy(plot_dom,error_ar_4,'--k','LineWidth',3);

h_legend = legend('10','5,5','2,2,2,2,2');%,'10,10');

%=========================
xlabel('Number of Coefficients', 'FontSize', 20) 
ylabel('Error in norm', 'FontSize', 20) 

set(gcf, 'Position', get(0,'ScreenSize'));
set(gca,'FontSize',30)
xlim([plot_dom(1),plot_dom(end)]);

if saveit
    saveas(gcf,name_it,'fig');
    saveas(gcf,name_it,'jpg');
    print('-depsc2',name_it);
end
