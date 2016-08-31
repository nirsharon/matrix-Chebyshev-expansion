function cjs=chebcoefs(func, n)
% Chebyshev fit: Given a function func, approximation such that 
%       f(x) \approx \sum_{k=0}^{n-1} c_{k}T_{k}(y) - c_{0}/2
% 
%
% Translated from Numerrical Recipes, Third edition, Section 5.8, pp. 236.

a =-1;
b =1;
if ~exist('n','var')
    n=50;
end

fks=zeros(n,1);
cjs=zeros(n,1);

bma=0.5*(b-a);
bpa=0.5*(b+a);

y=zeros(n,1);

for k=0:n-1 % We evaluate the function at the n points required by (5.8.7).
    y(k+1)=cos(pi*(k+0.5)/n); 
    fks(k+1)=func(y(k+1)*bma+bpa);
end
fac=2.0/n;
for j=0:n-1 % Now evaluate (5.8.7).
    sum=0.0;
    for k=0:n-1
        sum = sum + fks(k+1)*cos(pi*j*(k+0.5)/n);
    end
    cjs(j+1)=fac*sum;
end

end