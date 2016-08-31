function y= FilterFunc(x,r0,risetime,ctr)
%
% Construct a window function centered around ctr whose radius (measured
% betweetn the two values of x for which the windows gets the values 1/2)
% is r0, and whose raise time from 0.01 to 0.999 is (to a very good
% approximation) is 2*risetime.
%
y=0.5*(1-erf(2/risetime*(abs(x-ctr)-r0)));
end
