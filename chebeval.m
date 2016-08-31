function fv=chebeval(coef,A,m)
% Chebyshev evaluation: The Chebyshev polynomial
%       \sum_{k=0}^{m-1} c_{k}T_{k}(A) - c_{0}/2*I
%
% Based on Numerrical Recipes, Third edition, Section 5.8, pp. 237.
%

%if issparse(A)
    if max(eigs(A),1)>1
        error('The matrix does not fit the approximation');
    end
% else
%     if norm(A)>1
%         error('The matrix does not fit the approximation');
%     end
% end

if m>numel(coef) % Lowest order of apprixmation is one, which corrponds to the constant approximation
    error('Approximation order is too high for the precomputed C');
end

if m<1
    error('Approximation order must be greater than 1');
end

d  = zeros(size(A));
dd = zeros(size(A));
n  = size(A,1);
y2 = 2*A; 

for j=m-1:-1:1 % Clenshaw’s recurrence.
    sv = d;
    d  = y2*d-dd+coef(j+1)*eye(n);
    dd = sv;
end
fv= A*d-dd+0.5*coef(1)*eye(n);
end
