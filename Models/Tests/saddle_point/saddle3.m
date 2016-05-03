%% Extrema points in distribution
% Extrema points are given by the equation: 
% n_ratio = - (x - mu) / (x + mu) e^(8 x mu)
% Unable to solve it analytically so instead solve numerically for a sample
% of different (mu, n_ratio) points. For certain (mu, n_ratio) points x can 
% take on several values (1 value or 3 values), due to the non-linearity of
% the equation. The three values represent: max, min, and saddle points. 
% For each of these three we interpolate the value of x inbetween the 
% sampled values of x. Exports the extrepolation function of a
% (mu, n_ratio) point.


%% SOLVER

% Size of sample ? number of (mu, n_ratio) points to solve for x.
ss = 10000; % 4000;
% We run the numerically solver 3 times for each (mu, n_ratio) point,
% hopefully catching all the 3 values that x might take.
repeat = 3;

% Uniformly randomly draw (mu, n_ratio) points from their respective ranges
g_mu = rand(ss,1) * 1.5; % mu in [0, 1.5]
g_n_ratio = 1 + rand(ss,1); % n_ratio in [1, 2]

% Vector with the x solutions from the solver.
g_x = NaN(repeat*ss,1);

% The level of detail used to interpolate between the sample values of x.
% 100^2 = 10.000 points evenly spread in the parameter space of mu and
% n_ratio, ie. [0, 1.5]X[1, 2]
resolution = 100;   
% Create the mesh used to interpolate. 
[v_mu, v_n_ratio] = deal( 0:1.5/(resolution-1):1.5, 1:1/(resolution-1):2 );
[M_mu, M_n_ratio] = meshgrid(v_mu, v_n_ratio);

% Create the symbolic varible x used in the solver
syms x

% Create time stamp used to name exported files
timestamp = char(datetime('now'), 'yyyyMMdd_HHmmss');

h = waitbar(0,'Calculating...');
tic;
for i=1:ss
    mu = g_mu(i);
    n_ratio = g_n_ratio(i);
    
    % Solve 3 times for each (mu, n_ratio) point, using a uniformly random 
    % drawn starting point for finding solutions. Mean ideal points of 
    % subpopulations are -mu and mu respectively, so all extrema must lie 
    % within. Starting points are drawn from [-mu mu].
    sol = NaN(1,repeat);
    for j=1:repeat
        
        % Solve the equation
        sol_j = vpasolve( (x+mu)*n_ratio + (x-mu)*exp(8*x*mu), x, [-mu mu], 'random', true);
        
        % In case the solver cannot find solution set x to NaN.
        if ~isempty(sol_j)
            sol(j) = sol_j;
        else
            sol(j) = nan;
        end
    end
    
    % Save solutions
    g_x( repeat*(i-1)+1 : repeat*i ) = sol;
    
    waitbar(i/ss);
end
toc
close(h);

% Export solutions. 
output = table(repelem(1:ss,repeat)', repmat(1:repeat,1,ss)', repelem(g_mu,repeat), repelem(g_n_ratio,repeat), g_x, ...
    'VariableNames', {'sample' 'solution' 'mu' 'n_ratio' 'x'});
save(strcat('output_', timestamp,'.mat'), 'output');

% Warning if solver unable to find solution
num_nan = sum(isnan(g_x));
if num_nan>0
   warning('No solution in %s percent of all cases.',num2str(num_nan/ss*100))
end


%% IDENTIFY THE TYPE OF EXTREMA

% Convert solution vector to matrix (making it slightly easier to visually 
% debug and validate that conditions). Each row is a (mu, n_ratio) point.
g_x_matrix = reshape(g_x', repeat, ss)';

% Saddle points have non-negative x (due to the assymmetry in the range of
% n_ratio). Additionally saddle points fulfill the following condition:
% x < 1/2 sqrt(4 mu^2 - 1)
cond_saddle = logical((g_x_matrix < repmat(sqrt(4*g_mu.^2 - 1)/2, 1, 3)) .* (g_x_matrix >= 0));

% The remaning points are maximum points. The strictly positive x points
% correspond to the bimodal peak on the right side of zero. 
cond_right = logical(~cond_saddle .* (g_x_matrix > 0));

% The non-positive x points are either bimodal peak on the left side of zero,
% or peak of unimodal distribution - we lump these two together.
cond_left = logical(~cond_saddle .* (g_x_matrix <= 0));

%[cond_saddle cond_right cond_left]

% Convert back to vector
c_s = reshape(cond_saddle', repeat*ss, 1);
c_r = reshape(cond_right', repeat*ss, 1);
c_l = reshape(cond_left', repeat*ss, 1);

% NaN if condition is not fulfilled.
g_x_s = NaN(repeat*ss,1);
g_x_r = NaN(repeat*ss,1);
g_x_l = NaN(repeat*ss,1);
% Value if condition is fulfilled
g_x_s( c_s ) = g_x( c_s );
g_x_r( c_r ) = g_x( c_r );
g_x_l( c_l ) = g_x( c_l );

%[g_x_s g_x_r g_x_l]


%% COLOURS
% Creating colors for scatter plot.

% Category for respectively mu > 0, 0.5, 1
g_mu_g = (0.5 < g_mu)/2 + (1 < g_mu)/2;
% Color conditions for mu
g_mu_c1 = (0.5 >= g_mu);
g_mu_c3 = (1 < g_mu);
% Set color for mu
g_mu_c = repmat([23 177 43]/255, ss, 1); % Color of for 0.5 < mu <= 1
g_mu_c(g_mu_c1,:) = repmat([243 94 90]/255, sum(g_mu_c1), 1); % Color of for mu <= 0.5
g_mu_c(g_mu_c3,:) = repmat([80 134 255]/255, sum(g_mu_c3), 1);  % Color of for mu > 1

% Category for respectively n_ratio > 1, 4/3, 5/3
g_n_ratio_g = round((1 + (4/3 < g_n_ratio)/3 + (5/3 < g_n_ratio)/3)*100)/100;
% Color conditions for n_ratio
g_n_ratio_c1 = (4/3 >= g_n_ratio);
g_n_ratio_c3 = (5/3 < g_n_ratio);
% Set color for n_ratio
g_n_ratio_c = repmat([23 179 183]/255, ss, 1); % Color of for 4/3 < n_ratio <= 5/3
g_n_ratio_c(g_n_ratio_c1,:) = repmat([193 133 6]/255, sum(g_n_ratio_c1), 1); % Color of for n_ratio <= 4/3
g_n_ratio_c(g_n_ratio_c3,:) = repmat([252 65 193]/255, sum(g_n_ratio_c3), 1); % Color of for n_ratio > 5/3


%% MISC

% Extend vector of mu and n_ratio so they fit length of g_x
g_mu_rep = repelem(g_mu,repeat);
g_n_ratio_rep = repelem(g_n_ratio,repeat);
g_mu_g_rep = repelem(g_mu_g,repeat,1);
g_n_ratio_g_rep = repelem(g_n_ratio_g,repeat,1);


%% INTERPOLATE

% All points lump together
x_interpol = scatteredInterpolant(g_mu_rep, g_n_ratio_rep, g_x);
z = x_interpol(M_mu, M_n_ratio);

% Remove NaN
g_x_s_nan = ~isnan(g_x_s);
g_x_r_nan = ~isnan(g_x_r);
g_x_l_nan = ~isnan(g_x_l);

x_s_interpol = scatteredInterpolant(g_mu_rep(g_x_s_nan), g_n_ratio_rep(g_x_s_nan), g_x_s(g_x_s_nan),'linear','nearest');
x_s = x_s_interpol(M_mu, M_n_ratio);

x_r_interpol = scatteredInterpolant(g_mu_rep(g_x_r_nan), g_n_ratio_rep(g_x_r_nan), g_x_r(g_x_r_nan),'linear','nearest');
x_r = x_r_interpol(M_mu, M_n_ratio);

x_l_interpol = scatteredInterpolant(g_mu_rep(g_x_l_nan), g_n_ratio_rep(g_x_l_nan), g_x_l(g_x_l_nan),'linear','nearest');
x_l = x_l_interpol(M_mu, M_n_ratio);

% Export interpolation 
save(strcat('extrema_', timestamp, '_ss', num2str(ss),'.mat'), 'x_s_interpol', 'x_r_interpol', 'x_l_interpol');


%% PLOT


% All 3 solutions lumped together

figure(1); % (x, mu) diagram
gscatter( g_x , g_mu_rep, g_n_ratio_g_rep);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
% The break point between max and saddle point.
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'k');
% Left and right subpop ideal points
plot([0, -1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % left subpop. ideal point, x_l = -mu
plot([0, 1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % right subpop. ideal point,x_r = mu
% The limited of the break point between unimodal and bimodal (n_ratio = 1)
%plot([-1.5, 1.5], [0.5, 0.5], 'Color',[1,1,1]); %,[0.9,0.9,0.9]); 
%plot([-1.5, 1.5], [1, 1], 'Color',[0.9,0.9,0.9]); 
hold off;

figure(2); % (x, n_ratio) diagram
gscatter( g_x , g_n_ratio_rep, g_mu_g_rep);
%scatter( g_f4 , repelem(g_n_ratio,repeat), 8, repelem(g_mu,repeat), 'filled');
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');
%c = colorbar('southoutside');
%c.Label.String = '\mu';

figure(3);
meshc(v_mu, v_n_ratio, z);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, z);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;



% Saddle points -- solution

figure(1);
gscatter( g_x_s , g_mu_rep, g_n_ratio_g_rep);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
% The break point between max and saddle point.
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'k');
% Left and right subpop ideal points
plot([0, -1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % left subpop. ideal point, x_l = -mu
plot([0, 1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % right subpop. ideal point,x_r = mu
% The limited of the break point between unimodal and bimodal (n_ratio = 1)
%plot([-1.5, 1.5], [0.5, 0.5], 'Color',[1,1,1]); %,[0.9,0.9,0.9]); 
%plot([-1.5, 1.5], [1, 1], 'Color',[0.9,0.9,0.9]); 
hold off;

figure(2);
gscatter( g_x_s , g_n_ratio_rep, g_mu_g_rep);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');

figure(3);
%meshc(v_mu, v_n_ratio, x_s);
meshc(v_mu, v_n_ratio, x_s);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, x_s);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;




% Right peak of bimodal distribution

figure(1);
gscatter( g_x_r , g_mu_rep, g_n_ratio_g_rep);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
% The break point between max and saddle point.
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'k');
% Left and right subpop ideal points
plot([0, -1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % left subpop. ideal point, x_l = -mu
plot([0, 1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % right subpop. ideal point,x_r = mu
% The limited of the break point between unimodal and bimodal (n_ratio = 1)
%plot([-1.5, 1.5], [0.5, 0.5], 'Color',[1,1,1]); %,[0.9,0.9,0.9]); 
%plot([-1.5, 1.5], [1, 1], 'Color',[0.9,0.9,0.9]); 
hold off;

figure(2);
gscatter( g_x_r , g_n_ratio_rep, g_mu_g_rep);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');

figure(3);
meshc(v_mu, v_n_ratio, x_r);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, x_r);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;


% Left peak of bimodal distribution and peak of unimodal

figure(1);
gscatter( g_x_l , g_mu_rep, g_n_ratio_g_rep);
xlim([-1.5 1.5]); ylim([0 1.5]);
xlabel('x'); ylabel('\mu');
hold on;
% The break point between max and saddle point.
secdiffmu = (0.50:0.01:1.5)';
secdiffx = 1/2 * sqrt( 4 .* secdiffmu.^2 - 1 );
plot(secdiffx,secdiffmu, 'k');
% Left and right subpop ideal points
plot([0, -1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % left subpop. ideal point, x_l = -mu
plot([0, 1.5], [0, 1.5], 'Color',[0.9,0.9,0.9]); % right subpop. ideal point,x_r = mu
% The limited of the break point between unimodal and bimodal (n_ratio = 1)
%plot([-1.5, 1.5], [0.5, 0.5], 'Color',[1,1,1]); %,[0.9,0.9,0.9]); 
%plot([-1.5, 1.5], [1, 1], 'Color',[0.9,0.9,0.9]); 
hold off;

figure(2);
gscatter( g_x_l , g_n_ratio_rep, g_mu_g_rep);
xlim([-1.5 1.5]); ylim([1 2]);
xlabel('x'); ylabel('n_l/n_r');

figure(3);
meshc(v_mu, v_n_ratio, x_l);
xlabel('\mu'); ylabel('n_l/n_r');

figure(4);
imagesc(v_mu, v_n_ratio, x_l);
set(gca,'ydir', 'normal');
xlabel('\mu'); ylabel('n_l/n_r'); zlabel('x');
colorbar;





