
res = 1000;

%g_mu = rand(res,1) * 1.5; % mu in [0,1.5]
g_mu = normrnd(0.63, 0.1, res, 1); % mu in [0,1.5]
outside_range = logical( (g_mu < 0) + (1.5 > g_mu) );
g_mu(outside_range) = rand(sum(outside_range),1) * 1.5;

g_n_ratio = 1 + rand(res,1); % n_ratio in [1,2]

g_f = NaN(res,1);

h = waitbar(0,'Calculating...');
tic;
for r=1:res
    mu = g_mu(r);
    n_ratio = g_n_ratio(r);
    syms x
        
    % Start search at x=0
    sol = vpasolve( (x+mu)*n_ratio + (x-mu)*exp(8*x*mu), x, 0);
    
    % Random start within centered 90% of the range
    %sol = vpasolve( (x+mu)*n_ratio + (x-mu)*exp(8*x*mu), x, [-mu mu]*0.90, 'random', true);
    
    % In case the solver can't find a solution.
    if ~isempty(sol)
        g_f(r) = sol;
    else
        g_f(r) = nan;
    end
    
    waitbar(r/res);
end
toc
close(h);

output = table((1:res)', g_mu, g_n_ratio, g_f, ...
    'VariableNames', {'id' 'mu' 'n_ratio' 'x'});

num_nan = sum(isnan(g_f));
if num_nan>0
   warning('No solution in %s percent of cases.',num2str(num_nan/res*100))
end

%% PLOT

% Fixing colors in scatter plot.
g_mu_g = (0.5 < g_mu)/2 + (1 < g_mu)/2;
g_mu_c = repmat([23 177 43]/255, res, 1);
g_mu_c1 = (0.5 >= g_mu);
g_mu_c(g_mu_c1,:) = repmat([243 94 90]/255, sum(g_mu_c1), 1);
g_mu_c3 = (1 < g_mu);
g_mu_c(g_mu_c3,:) = repmat([80 134 255]/255, sum(g_mu_c3), 1);

g_n_ratio_g = 1 + (4/3 < g_n_ratio)/3 + (5/3 < g_n_ratio)/3;
g_n_ratio_c = repmat([23 179 183]/255, res, 1);
g_n_ratio_c1 = (4/3 >= g_n_ratio);
g_n_ratio_c(g_n_ratio_c1,:) = repmat([193 133 6]/255, sum(g_n_ratio_c1), 1);
g_n_ratio_c3 = (5/3 < g_n_ratio);
g_n_ratio_c(g_n_ratio_c3,:) = repmat([252 65 193]/255, sum(g_n_ratio_c3), 1);

% Mesh

s = 100;
[v_mu, v_n_ratio] = deal( 0:1.5/(s-1):1.5, 1:1/(s-1):2 );
[M_mu, M_n_ratio] = meshgrid(v_mu, v_n_ratio);


% First

figure(1);
gscatter( g_f , g_mu, g_n_ratio_g);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
plot([-1.5, 1.5], [0.63, 0.63], 'k--');
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'm:');
hold off;

figure(2);
gscatter( g_f , g_n_ratio, g_mu_g);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');

x_interpol1 = scatteredInterpolant(g_mu, g_n_ratio, g_f);
z1 = x_interpol1(M_mu, M_n_ratio);

figure(3);
meshc(v_mu, v_n_ratio, z1);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, z1);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');



% Tests

outliers1 = (g_f > 0.4);
[g_mu(outliers1) g_n_ratio(outliers1) g_f(outliers1)]
mean([g_mu(outliers1) g_n_ratio(outliers1) g_f(outliers1)])
std([g_mu(outliers1) g_n_ratio(outliers1) g_f(outliers1)])



