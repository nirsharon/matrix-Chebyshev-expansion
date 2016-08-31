% script name: "non_smooth_f1"
% Testing the first function: sign(x)*x^2

saveit = 1;
name_it = 'non_smooth_f1';

k        = 10;  % matrix size
coef_num = 2000; %50;

% the scalar function
f = @(x) sign(x)*x^2;
coefs = chebcoefs(f,coef_num);
coefs(1:2:coef_num) = 0;

% random symmetric matrix
A = rand(k);
A = A+A'   ;
A = A/norm(A,'fro');   % verify spectrum in [-1,1]

% matrix function by spectral decmposition, our ground truth
[V, D] = eig(A);
f_A    =  V*(sign(D)*D^2)*V';  
cn     = cond(A)

% initialize arrays
number_of_plots = 20;
plot_dom        = floor(linspace(1,coef_num,number_of_plots));
error_ar        = zeros(1,number_of_plots);

% main loop
for l=1:number_of_plots
    Cheb = chebeval(coefs,A,plot_dom(l));
    error_ar(l) = norm(Cheb-f_A)/cn;
end;



%========== plot ==========
figure;
semilogy(2:2:coef_num,abs(coefs(2:2:coef_num)),'--k','LineWidth',4);
hold on
semilogy(plot_dom,error_ar,'LineWidth',4);
% theoretical bound
offset = plot_dom(1)^2*error_ar(1);
semilogy(plot_dom,offset*plot_dom.^(-2),'r','LineWidth',4.5);
xlabel('n');

%set(gca,'LineWidth',1.5);
set(gcf, 'Position', get(0,'ScreenSize'));
set(gca,'FontSize',30)
str1 = '\alpha_n[f] (absolute value)';
str2 = 'Truncation error (Frobenius norm)';
str3 = 'theoretical bound of n^{-2}';


h_legend = legend(str1,str2,str3);

%=========================
% y1 = 1.1*min(error_ar);
% y2 = 1.1*max(error_ar);
% ylim([y1 y2]);

if saveit
    saveas(gcf,name_it,'fig');
    saveas(gcf,name_it,'jpg');
    print('-depsc2',name_it);
end
