function [p, z, za, zb] = corr_rtest(ra, rb, na, nb)
%[p, z, za, zb] = corr_rtest(ra, rb, na, nb)
%
% comparison of two correlation coefficients
% 
% input: ra, rb: r-values of correlation coefficient
%        na, nb: the number of datas
%
% ouput: p(1): pvalue of one-tailed test
%        p(2): pvalue of two-tailed test
%        za, zb: z transformed value of ra, rb
%{
e.g.
x = 0:.01:3;

noise1 = 0.12*rand(size(x));%randn
noise2 = 0.01*rand(size(x));

y1 = 1.2*x + noise1+0;  %y1 = 4*x + noise1 + 5;
y2 = 1.2*x + noise2+0;  %y2 = 2*x + noise2 + 3;

figure(1); clf; hold on;
plot(x, y1, 'ob');
plot(x, y2, 'or');

r1 = corr(x', y1')
r2 = corr(x', y2')
[p, z] = corr_rtest(r1, r2, len(x), len(x))

%}

%Ryosuke F Takeuchi 2017/02/02 - 

	za  = r2z(ra);
	zb  = r2z(rb);
	szab = sqrt(1/(na-3) + 1/(nb-3));
	z = abs(za-zb)/szab;
	p(1) = 1-normcdf(z, 0, 1);
	p(2) = 2*p(1);	
end

function zval = r2z(rval)
  zval = 0.5*log((1+rval)/(1-rval));
end