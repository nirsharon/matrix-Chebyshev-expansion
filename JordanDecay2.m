% script name: "JordanDecay2"
% Compare the truncation error of three Jordan blocks of different sizes


% set "saveit" to 1 to save the results into "nameit.fig/jpf/eps" files
saveit = 0;
nameit = 'JordanBlocks';

deg  = 3;
func = @(x) abs(x).^(deg+.5); 
N    = 1000;
coef = chebcoefs(func,N);

% % the decay rate of coefficients
% figure; semilogy(1:2:N,abs(coef(1:2:N)));
% hold on; semilogy(1:2:N,(1:2:N).^(-deg-1),'r');

lambda     = .7;

number_of_truncations = 16;

N_sample   = round(linspace(10,N,number_of_truncations));
m_sizes    = [2,3,4];
error_rate = zeros(number_of_truncations,numel(m_sizes));
ErrorMat   = cell(numel(m_sizes),1);

for m_iter =1:numel(m_sizes)
    % construct the Jordan block of order m
    m = m_sizes(m_iter);
    J = lambda*eye(m)+triu(ones(m),1)-triu(ones(m),2);
    
    % our ground truth, we set lambda to be positive
    GT = sqrtm(J)*(J)^deg; 
    for j=1:number_of_truncations
        F_J = chebeval(coef,J,N_sample(j));
        error_rate(j,m_iter)= norm(F_J - GT,'fro')/cond(J)+eps;
        % The last error matrix
        if j==number_of_truncations
            ErrorMat{m_iter,1} = abs(F_J - GT);
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

str1 = ['m = ',num2str(m_sizes(plotind(3)))];
str2 = ['m = ',num2str(m_sizes(plotind(2)))];
str3 = ['m = ',num2str(m_sizes(plotind(1)))];
h_legend = legend(str1,str2,str3);

xlabel('Expansion length');
ylabel('Error (Frobenius norm)');
set(gca,'FontSize',26)
set(gcf, 'Position', get(0,'ScreenSize'));

%set(h_legend,'FontSize',27);

if saveit
    saveas(gcf,nameit,'fig');
    saveas(gcf,nameit,'jpg');
    print('-depsc2',nameit);
end
