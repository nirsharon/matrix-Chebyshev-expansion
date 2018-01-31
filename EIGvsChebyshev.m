% script name: "EIGvsChebyshev"

% global constants
halfs_len = 20;         % eigenspace size
center    = 0.5;        % eigenvalue value

% comparison setting
matrixSizes = round(2*logspace(3,4,10));
iterNum     = numel(matrixSizes);
Cheb_time   = zeros(iterNum,1);
eig_time    = zeros(iterNum,1);

% main loop
for count=1:iterNum
    disp(['Iteration number ', num2str(count), ' out of ',num2str(iterNum)]);
    n = matrixSizes(count);
    
    % construct the spectrum
    ones_len  = floor(n/2);
    zeros_len = n-ones_len-halfs_len;
    sigma  = [ones(ones_len,1) ; center*ones(halfs_len,1) ; zeros(zeros_len,1)];
    
    % the matrix-vector evaluation
    %AxFunc =  @(x) idctt(bsxfun(@times,sigma,dctt(x)));
    
    [B, ~] = qr(rand(n));
    A = B*diag(sigma)*B';
    AxFunc =  @(x) A*x;
    
    
    % the filter
    filterFunc = @(x) FilterFunc(x,.25,.15,center);
    cc         = chebcoefs(filterFunc,7000);
    m     = 10;
    coef  = cc(1:m);
    
    % comparison
    % EIGS
    opts.v0 = x;
    invAfun = @(x) (A-center*eye(n))\x;
    opts.tol = 1e-10;
    tic
    [ueig, eig_res, flagg] = eigs(A,halfs_len,center); %,opts);
    t = toc;
    eig_time(count) = t;    
    norm(AxFunc(ueig(:,1:halfs_len))-center*ueig(:,1:halfs_len))
    %norm(AxFunc(ueig(:,1:halfs_len))-ueig(:,1:halfs_len)*eig_res)
    if flagg
        disp('No convergence in EIGS');
    end
    
    x = rand(n,halfs_len);   % arbitrary initial guess
    warning off;
    precision = 1e-10;
    diff    = precision+1;
    maxIter = 20;
    iter    = 1;
    %Chebyshev
    tic
    res = chebeval_vector_AxFunc(coef, AxFunc ,m , x);
    % Repeating for higher accuracy
    while (diff>precision)&&(iter<maxIter)
        resn = chebeval_vector_AxFunc(coef, AxFunc ,m , res);
        diff=norm(AxFunc(resn)-center*resn,'fro')/norm(resn(:));
        %diff = norm(resn-res)/norm(res);
        res = resn;
        iter = iter +1;
    end
    t = toc;
    Cheb_time(count) = t;
    norm(AxFunc(res(:,1:halfs_len))-center*res(:,1:halfs_len))

    warning on;
end

% results
res = [Cheb_time , eig_time]
plot(matrixSizes,res,'LineWidth',4);
h = legend('Chebyshev','EIGS','Location','NorthWest'); %'northoutside','Orientation','horizontal')
grid on
set(gca,'FontSize',27);
xlabel('Matrix size');
ylabel('Time in seconds');

% % validate result only for small enough arrays
% if n<5000
%     [u,s,v] = svd(res);
%     bar(diag(s))
%     norm(AxFunc(u(:,1:halfs_len))-center*u(:,1:halfs_len))
%
%     if size(res,2)==size(res,1)
%         sigma_tilde = [zeros(ones_len,1) ; zeros(zeros_len,1); center*ones(halfs_len,1)];
%         norm(sort(sigma_tilde)-sort(real(eig(res))))
%     end
% end
%
%
%
