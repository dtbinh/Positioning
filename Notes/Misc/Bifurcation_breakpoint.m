% What is the value of mu for which the distribution changes from
% unimodal to bimodal for a given value of n_ratio.
n = 2
syms x
vpasolve( (sqrt(4*x^2-1)-2*x)/(sqrt(4*x^2-1)+2*x) * exp(4*sqrt(4*x^2-1)*x) + n == 0, x, 'random', true) 

% n = 2 gives mu = 0.65687728299045653134597517069296
% n = 1 gives mu = 0.50000000000000000000000000000775




% Extrema points 
extrema(0.5, 1)
extrema(1, 1)
extrema(1.5, 1)
extrema(0.5, 2)
extrema(1, 2)
extrema(2, 2)