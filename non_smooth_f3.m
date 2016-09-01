% script name: "non_smooth_f3"
% Testing the third function: 1/x^2+.25

saveit = 1;
name_it      = 'non_smooth_f3';
name_it_spec = 'spectrum_non_smooth_f3';


k        = 10;  % matrix size
coef_num = 73;

% the scalar function
f = @(x) (x.^2+.25).^(-1);
coefs = chebcoefs(f,coef_num);
coefs(2:2:coef_num) = 0;

% The matrix
A      = rand(k);
[q, ~] = qr(A);
spect = 2*rand(1,k)-1;
spect = spect/max(abs(spect));
spect = sort(spect);
A = q'*diag(spect)*q;
cn = cond(A)

% matrix function by applying inverse
f_A =  inv(A^2+.25*eye(k));

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
semilogy(1:2:coef_num,abs(coefs(1:2:coef_num)),'--k','LineWidth',4);
hold on
semilogy(plot_dom,error_ar,'LineWidth',4);
%xlabel('n');

set(gcf, 'Position', get(0,'ScreenSize'));
set(gca,'FontSize',26)
str1 = '\alpha_n[f] (absolute value)';
str2 = 'Truncation error (Frobenius)';
h_legend = legend(str1,str2);
xlim([plot_dom(1),plot_dom(end)]);


if saveit
    saveas(gcf,name_it,'fig');
    saveas(gcf,name_it,'jpg');
    print('-depsc2',name_it);
end


% spectrum plot
figure;
bar(spect);
hold on;
plot(0:11,ones(12,1)*(-.5),'r');
plot(0:11,ones(12,1)*(.5),'r');

y1 = min(spect);
y2 = max(spect);
ylim([y1-.1 y2+.1]);
xlim([0,11]);


%xlabel('The spectrum of the matrix','FontSize', 26)

set(gca,'LineWidth',1.5);
set(gcf, 'Position', get(0,'ScreenSize'));
set(gca,'FontSize',26)

if saveit
    saveas(gcf,name_it_spec,'fig');
    saveas(gcf,name_it_spec,'jpg');
    print('-depsc2',name_it_spec);
end

