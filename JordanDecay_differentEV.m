% script name: "JordanDecay_differentEV"
% Compare the truncation error of three Jordan blocks of different
% eigenvalues


% set "saveit" to 1 to save the results into "nameit.fig/jpf/eps" files
saveit = 1;
nameit = 'JordanBlocks_diffEV';

deg  = 3;
m    = 3;   % convergence guaranteed for m<=deg with lambda<1
func = @(x) abs(x).^(deg+.5); 
N    = 2000;
coef = chebcoefs(func,N);

% % the decay rate of coefficients
% figure; semilogy(1:2:N,abs(coef(1:2:N)));
% hold on; semilogy(1:2:N,(1:2:N).^(-deg-1),'r');

lambda_vals           = [.4, .7, 1];
number_of_truncations = 16;

N_sample   = round(linspace(10,N,number_of_truncations));
error_rate = zeros(number_of_truncations,numel(lambda_vals));
ErrorMat   = cell(numel(lambda_vals),1);

for l_iter =1:numel(lambda_vals)
    % construct the Jordan block of order m
    lambda = lambda_vals(l_iter);
    J      = lambda*eye(m)+triu(ones(m),1)-triu(ones(m),2);
    
    % our ground truth, we set lambda to be positive
    GT = sqrtm(J)*(J)^deg; 
    for j=1:number_of_truncations
        F_J = chebeval(coef,J,N_sample(j));
        error_rate(j,l_iter)= norm(F_J - GT,'fro')/cond(J)+eps;
        % The last error matrix
        if j==number_of_truncations
            ErrorMat{l_iter,1} = abs(F_J - GT);
        end
    end
end

plotind = [1,2,3];
% plotting

figure; 
semilogy(N_sample,error_rate(:,(plotind(3))),'b','LineWidth',4);
hold on;
semilogy(N_sample,error_rate(:,(plotind(2))),':k','LineWidth',3.5);
semilogy(N_sample,error_rate(:,(plotind(1))),'--r','LineWidth',4.5);

str1 = ['\lambda = ',num2str(lambda_vals(plotind(3)))];
str2 = ['\lambda = ',num2str(lambda_vals(plotind(2)))];
str3 = ['\lambda = ',num2str(lambda_vals(plotind(1)))];
h_legend = legend(str1,str2,str3,'Location','best');
 


xlabel('Expansion length');
ylabel('Error (Frobenius norm)');
set(gca,'FontSize',30)
set(gcf, 'Position', get(0,'ScreenSize'));

set(h_legend,'FontSize',28);

if saveit
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end
