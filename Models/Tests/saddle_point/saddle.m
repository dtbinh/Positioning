
res = 4000;
%[g_mu, g_n_ratio] = deal( 0 : 1.5/(res-1) : 1.5, 1 : 1/(res-1) : 2 );

g_mu = rand(res,1) * 1.5; % mu in [0,1.5]
g_n_ratio = 1 + rand(res,1); % n_ratio in [1,2]
%g_f = NaN(res,1);
repeat = 3;
g_f4 = NaN(repeat*res,1);
h = waitbar(0,'Calculating...');
tic;
for r=1:res
    mu = g_mu(r);
    n_ratio = g_n_ratio(r);
    %syms f(x)
    %f(x) = (x+mu)*n_ratio + (x-mu)*exp(8*x*mu);
    syms x
    
    % Preform 4 repetions with random start. Select the x-point closest to
    % zero.
    %sol = NaN(1,4);
    %for i=1:4
    %    sol(i) = vpasolve( (x+mu)*n_ratio + (x-mu)*exp(8*x*mu), x ,[-mu mu], 'random', true);
    %end
    %[~, min_i] = min(abs(sol));
    %g_f(r) = sol(min_i);
    
    % Random start.
    %g_f(r) = vpasolve( (x+mu)*n_ratio + (x-mu)*exp(8*x*mu), x ,[-mu mu],'random',true);
    
    % Preform 4 repetions with random start. Save all.
    sol = NaN(1,repeat);
    for i=1:repeat
        sol_i = vpasolve( (x+mu)*n_ratio + (x-mu)*exp(8*x*mu), x, [-mu mu], 'random', true);
        %sol_i = vpasolve( (x+mu)*n_ratio + (x-mu)*exp(8*x*mu), x, 0);
        
        % In case the solver can't find a solution.
        if ~isempty(sol_i)
            sol(i) = sol_i;
        else
            sol(i) = nan;
        end
    end
    g_f4((repeat*(r-1)+1):repeat*r) = sol;
    
    waitbar(r/res);
end
toc
close(h);

output = table(repelem(1:res,repeat)', repmat(1:repeat,1,res)', repelem(g_mu,repeat), repelem(g_n_ratio,repeat), g_f4, ...
    'VariableNames', {'id' 'solution' 'mu' 'n_ratio' 'x'});

num_nan = sum(isnan(g_f4));
if num_nan>0
   warning('No solution in %s percent of cases.',num2str(num_nan/res*100))
end

%% ALTERNATIVE SOLUTIONS:

% Only the first solution the solver finds.
g_f1 = g_f4(1:repeat:res*repeat);

% Only plot the solution closest to 0.
g_f4_matrix = reshape(g_f4', repeat, res)';
[~, min_i] = min(abs(g_f4_matrix), [], 2);
g_fm = g_f4_matrix( sub2ind(size(g_f4_matrix),(1:res)',min_i) );

% Plot the lowest x solution.
g_flow = min(g_f4_matrix, [], 2);

% Plot saddle points
% Using x < 1/2 sqrt(4 mu^2 - 1) and x >= 0
[condition, condition_i] = max((g_f4_matrix < repmat(sqrt(4*g_mu.^2 - 1)/2, 1, 3)) .* (g_f4_matrix >= 0), [], 2);
condition = logical(condition);
condition_j = (1:res)';
condition_ind = sub2ind(size(g_f4_matrix), condition_j((condition)), condition_i((condition)) );
g_fs = NaN(res,1); % Those that don't furfill condition is set to nan
g_fs((condition)) = g_f4_matrix( condition_ind );

% negation of condition
%g_fNOTs((~condition)) = g_f4_matrix( setdiff((1:res*repeat)', condition_ind) );


%% PLOT

% Fixing colors in scatter plot.
g_mu_g = (0.5 < g_mu)/2 + (1 < g_mu)/2;
g_mu_c = repmat([23 177 43]/255, res, 1);
g_mu_c1 = (0.5 >= g_mu);
g_mu_c(g_mu_c1,:) = repmat([243 94 90]/255, sum(g_mu_c1), 1);
g_mu_c3 = (1 < g_mu);
g_mu_c(g_mu_c3,:) = repmat([80 134 255]/255, sum(g_mu_c3), 1);

g_n_ratio_g = round((1 + (4/3 < g_n_ratio)/3 + (5/3 < g_n_ratio)/3)*100)/100;
g_n_ratio_c = repmat([23 179 183]/255, res, 1);
g_n_ratio_c1 = (4/3 >= g_n_ratio);
g_n_ratio_c(g_n_ratio_c1,:) = repmat([193 133 6]/255, sum(g_n_ratio_c1), 1);
g_n_ratio_c3 = (5/3 < g_n_ratio);
g_n_ratio_c(g_n_ratio_c3,:) = repmat([252 65 193]/255, sum(g_n_ratio_c3), 1);

%figure(1);
%scatter( g_f , g_mu, 15, g_n_ratio_c, 'filled');
%figure(2);
%scatter( g_f , g_n_ratio, 15, g_mu_c, 'filled');


% Mesh

s = 100;
[v_mu, v_n_ratio] = deal( 0:1.5/(s-1):1.5, 1:1/(s-1):2 );
[M_mu, M_n_ratio] = meshgrid(v_mu, v_n_ratio);


% All four

figure(1);
gscatter( g_f4 , repelem(g_mu,repeat), repelem(g_n_ratio_g,repeat,1));
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
plot([-1.5, 1.5], [0.5, 0.5], 'k--');
plot([-1.5, 1.5], [1, 1], 'k--'); % 
plot([0, -1.5], [0, 1.5], 'k--'); % x = -mu
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'm:');
hold off;

figure(2);
gscatter( g_f4 , repelem(g_n_ratio,repeat), repelem(g_mu_g,repeat,1));
%scatter( g_f4 , repelem(g_n_ratio,repeat), 8, repelem(g_mu,repeat), 'filled');
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');
%c = colorbar('southoutside');
%c.Label.String = '\mu';

x_interpol = scatteredInterpolant(repelem(g_mu,repeat), repelem(g_n_ratio,repeat), g_f4);
z = x_interpol(M_mu, M_n_ratio);

figure(3);
meshc(v_mu, v_n_ratio, z);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, z);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;

% First

figure(1);
gscatter( g_f1 , g_mu, g_n_ratio_g);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
plot([-1.5, 1.5], [0.5, 0.5], 'k--');
plot([-1.5, 1.5], [1, 1], 'k--');
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'm:');
hold off;

figure(2);
gscatter( g_f1 , g_n_ratio, g_mu_g);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');

x_interpol1 = scatteredInterpolant(g_mu, g_n_ratio, g_f1);
z1 = x_interpol1(M_mu, M_n_ratio);

figure(3);
meshc(v_mu, v_n_ratio, z1);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, z1);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;

% Min / closest to zero

figure(1);
gscatter( g_fm , g_mu, g_n_ratio_g);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
plot([-1.5, 1.5], [0.5, 0.5], 'k--');
plot([-1.5, 1.5], [1, 1], 'k--');
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'm:');
hold off;

figure(2);
gscatter( g_fm , g_n_ratio, g_mu_g);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');

x_interpolm = scatteredInterpolant(g_mu, g_n_ratio, g_fm);
zm = x_interpolm(M_mu, M_n_ratio);

figure(3);
meshc(v_mu, v_n_ratio, zm);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, zm);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;


% Lowest solution

figure(1);
gscatter( g_flow , g_mu, g_n_ratio_g);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
plot([-1.5, 1.5], [0.5, 0.5], 'k--');
plot([1, 1], [0, 1.5], 'k--');
plot([-1, -1], [0, 1.5], 'k--');
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'm:');
hold off;

figure(2);
gscatter( g_flow , g_n_ratio, g_mu_g);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');

x_interpollow = scatteredInterpolant(g_mu, g_n_ratio, g_flow);
zlow = x_interpollow(M_mu, M_n_ratio);

figure(3);
meshc(v_mu, v_n_ratio, zlow);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, zlow);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;



% Saddle points -- solution

figure(1);
gscatter( g_fs , g_mu, g_n_ratio_g);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
plot([-1.5, 1.5], [0.5, 0.5], 'k--');
plot([1, 1], [0, 1.5], 'k--');
plot([-1, -1], [0, 1.5], 'k--');
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'm:');
hold off;

figure(2);
gscatter( g_fs , g_n_ratio, g_mu_g);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');


g_fs_nan = ~isnan(g_fs);
x_s_interpol = scatteredInterpolant(g_mu(g_fs_nan), g_n_ratio(g_fs_nan), g_fs(g_fs_nan),'linear','nearest')
x_s = x_s_interpol(M_mu, M_n_ratio);

%M_zs = griddata(g_mu(g_fs_nan), g_n_ratio(g_fs_nan), g_fs(g_fs_nan), M_mu, M_n_ratio);

figure(3);
%meshc(v_mu, v_n_ratio, zs);
mesh(v_mu, v_n_ratio, x_s);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, x_s);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;

%figure(5);
%mesh(M_mu, M_n_ratio, M_zs);
%xlim([0 1.5]); ylim([1 2]); %zlim([-1.5 1.5]);
%xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');









