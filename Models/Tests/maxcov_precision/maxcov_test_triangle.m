
clearvars;

%rng('default');

pref.repetitions = 200;

%% Boundary box

% Boundary box
%bbox = [x([1 end end 1 1]); y([end end 1 1 end])]';
bbox = [-5 5 5 -5 -5; 5 5 -5 -5 5]';

% Boundary box line segments
bvx = [ bbox(1:end-1,1)'; circshift( bbox(1:end-1,1)', [0,-1] ) ];
bvy = [ bbox(1:end-1,2)'; circshift( bbox(1:end-1,2)', [0,-1] ) ];

%% Population
pref.boundary = 12; % Number of standard deviations
pref.resolution = 50; % Length of the square (even number to include (0,0))

sd = 0.5; % Standard deviation of each subpopulation
b = pref.boundary/2;

% Grid / Square
[x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
[X,Y] = meshgrid(x,y);


%% Repetitions
xy_test = NaN(12, 2, pref.repetitions);
%correct = NaN(pref.repetitions,1);
error_share = NaN(pref.repetitions,1);

for rep=1:pref.repetitions
    
    pref.N = 10;
    pref.mu = 0; % Mean of subpopulation
    pref.n_ratio = 1; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
    
    %% Setup
    % Points - Random firm locations
    [x0, y0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*3 );
    xy_test(1:pref.N,:,rep) = [x0 y0];
    
    
    %% Population
    
    % Subpopulation right
    mu_r = [pref.mu 0]; % Subpopulation mean only deviate on x-axis
    sigma_r = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
    F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r); % Subpopulation pdf evaluated at each point in grid/square
    F_r = reshape(F_r,length(y),length(x)); % Formated as grid

    % Subpopulation left
    mu_l = [-pref.mu 0]; % Subpopulation mean only deviate on x-axis
    sigma_l = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
    F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l); % Subpopulation pdf evaluated at each point in grid/square
    F_l = reshape(F_l,length(y),length(x)); % Formated as grid

    % Total population 
    % The left subpopulation share of total population
    weight = pref.n_ratio/(1+pref.n_ratio);
    % Population probability density function (PDF). Left subpopulation is 
    % left unchanged, while the right subpopulation will scale according to 
    % n_ratio. Deviding by 2^2 has no effect on results, but scales PDF so 
    % its comparable to the baseline distribution with mean (0,0) and 
    % std. dev. (1,1). The subpopulations halves the std dev., thus their 
    % variance will be one-fourth.
    F = (F_l + F_r*1/pref.n_ratio)/4;
    
    
    %% Maxcov
    [xy_centroid_tri, xy_new_idx] = maxcov_delaunay_test(xy_test(1:pref.N,:,rep), [X(:) Y(:)], F);

    %% Market share
    triangles = size(xy_centroid_tri,1);
    shares = NaN(triangles, 1);
    for tri=1:triangles
        xy_tri = [xy_test(1:pref.N,:,rep); xy_centroid_tri(tri,:)];
        [market_i, ~] = marketshare4(xy_tri, [X(:) Y(:)]);
        
        % Index for each customer of the firm.
        idx = find(market_i==pref.N+1);
        
        shares(tri) = sum( F(idx) )/sum(F(:));
    end
    [share_max, share_max_idx] = max(shares);
    
    %correct(rep) = (xy_new_idx == share_max_idx);
    error_share(rep) = share_max-shares(xy_new_idx);
    error_distance(rep) = pdist([xy_centroid_tri(xy_new_idx,:); xy_centroid_tri(share_max_idx,:)], 'euclidean');
end


num_triangles = pref.N*2+2;
random_corrent_pct = 1/num_triangles;

correct = (error_share == 0);
corrent_pct = sum(correct)./pref.repetitions;
[max_error_share, max_error_share_rep] = max(error_share);
max_error_distance = error_distance(max_error_share_rep);

mean_error_share = mean(error_share(~correct));
mean_error_distance = mean(error_distance(~correct));

summary = table(random_corrent_pct, corrent_pct, max_error_share, mean_error_share, max_error_distance, mean_error_distance)





xy = xy_test(1:pref.N,:,max_error_share_rep);

figure(10);
%histogram(error_share);
histfit(error_share(~correct), [], 'exponential');
title('distribution of share errors');
xlabel('error in share (pct.)');

figure(11);
histfit(error_distance(~correct), [], 'Lognormal');
title('distribution of distance errors');
xlabel('error in distance (std. dev.)');