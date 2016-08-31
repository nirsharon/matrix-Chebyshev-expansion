% script name: "HugeMatrixExample"
%
%
% This script runs the example of eigenspace recovery for FFT matrix

% global constants
halfs_len = 20;         % eigenspace size
center    = 0.5;        % eigenvalue value

% comparison setting
matrixSizes = [5*1e4, 1e5, 5*1e5, 1e6]; %round(2*logspace(3,4,10));
iterNum     = numel(matrixSizes);
Cheb_time   = zeros(iterNum,1);
diff        = zeros(iterNum,1);

fid=fopen('rest.txt','w');

% main loop
for count=1:iterNum
    disp(['Iteration number ', num2str(count), ' out of ',num2str(iterNum)]);
    n = matrixSizes(count);
    
    % construct the spectrum
    ones_len  = floor(n/2);
    zeros_len = n-ones_len-halfs_len;
    sigma  = [ones(ones_len,1) ; center*ones(halfs_len,1) ; zeros(zeros_len,1)];
    
    % the matrix-vector evaluation
    AxFunc =  @(x) idctt(bsxfun(@times,sigma,dctt(x)));
    
%     [B, ~] = qr(rand(n));
%     A = B*diag(sigma)*B';
%     AxFunc =  @(x) A*x;
%     
    
    % the filter
    filterFunc = @(x) FilterFunc(x,.25,.15,center);
    cc         = chebcoefs(filterFunc,7000);
    m     = 10;
    coef  = cc(1:m);   
    
    x = rand(n,halfs_len);   % arbitrary initial guess
    precision = 1e-10;
    cdiff    = precision+1;
    maxIter = 20;
    iter    = 1;
    %Chebyshev
    tic
    res = chebeval_vector_AxFunc(coef, AxFunc ,m , x);
    % Repeating for higher accuracy
    while (cdiff>precision)&&(iter<maxIter)
        resn = chebeval_vector_AxFunc(coef, AxFunc ,m , res);
        cdiff=norm(AxFunc(resn)-center*resn,'fro')/norm(resn(:));
        %diff = norm(resn-res)/norm(res);
        res = resn;
        iter = iter +1;
    end
    t = toc;
    Cheb_time(count) = t
    diff(count) =norm(AxFunc(res(:,1:halfs_len))-center*res(:,1:halfs_len))
    iter
    
    fprintf(fid,'%1.1e & %5.3e & %7.3f  \\\\\n',n,diff(count),t);

end
fclose(fid);
fid=fopen('horiz_table.txt','w');
    fprintf(fid,'size          & %1.1e & %1.1e & %1.1e & %1.1e  \\\\\n',matrixSizes);
    fprintf(fid,'Timing (sec.) & %7.3f & %7.3f & %7.3f & %7.3f  \\\\\n',Cheb_time(:)  );
    fprintf(fid,'Precision     & %5.3e & %5.3e & %5.3e & %5.3e   \\\\\n',diff (:));
fclose(fid);
