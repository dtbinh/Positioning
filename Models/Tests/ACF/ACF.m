
% Autocorrelation function
%x = [-1 1 -1 1];
%x = [1 1.1 0.9 1.2];
x = cos(rand(1, 4)*2*pi-pi)

%x_0 = conj(x); x_1 = conj(x); x_2 = conj(x); x_3 = conj(x);
%x_1(1) = [];
%x_2(1:2) = [];
%x_3(1:3) = [];

% Lags:
%p_0 = sum((x-mean(x)) .* (x_0-mean(x))) / sum((x-mean(x)).^2);
%p_1 = sum((x(1:end-1)-mean(x)) .* (x_1-mean(x))) / sum(x.^2);
%p_2 = sum(x(1:end-2).*x_2)/sum(x.^2);
%p_3 = sum(x(1:end-3).*x_3)/sum(x.^2);

%[p_0 p_1 p_2 p_3]
[lags, ~, ~] = autocorr(x)

p2_0 = p_k(x, 0);
p2_1 = p_k(x, 1);
p2_2 = p_k(x, 2);
p2_3 = p_k(x, 3);

lags2 = [p2_0 p2_1 p2_2 p2_3]


% confidence interval is plus-minus z/sqrt(N). At 85% confidence z = 1.44.
z = 1.022;
ci = z/sqrt(length(x));
ci2 = -1/length(x) - z/sqrt(length(x));
ci_anderson1m = ( -1 - z * sqrt(length(x)-1-1) ) / (length(x)-1);
ci_anderson1p = ( -1 + z * sqrt(length(x)-1-1) ) / (length(x)-1);
%ci4 = exp(2*z)/(sqrt(length(x)-3)-1) / exp(2*z)/(sqrt(length(x)-3)+1);

%se_p1 = sqrt( 1/length(x) * (1 - 2 * sum([p2_0 + p2_1].^2)) );

confidence = [ci ci2 ci_anderson1m ci_anderson1p]

figure(1);
plot(x, 'b*');
ylim([-1 1]);

figure(3);
clf reset; % Reset figure.
stem(0:3, lags, 'MarkerEdgeColor', 'b', 'Color', 'b');
hold on;
stem(0:3, lags2, 'MarkerEdgeColor', 'g', 'Color', 'g');
scatter([1 1], [ci_anderson1m ci_anderson1p], '.');
hold off;
ylim([-1.2 1.2]);
