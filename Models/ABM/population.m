%% Population
% version 0.01
% Jonas K. Sekamane
%
% Inspired in part by: 
%   Lever and Sergenti (2011)

clearvars;

%% Preferences
rng('default'); % Seed such that the randomly generated results are repeatable
pref.boundary = 10; % Number of standard deviations
pref.resolution = 20; % Length of the square (even number to include (0,0))
pref.N = 5; % Number of firms
pref.mu = 1.5; % Mean of subpopulation
pref.n_ratio = 2; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
pref.iterations = 24;

%% Setup
unit = pref.resolution/pref.boundary; % 1 standard deviation
firms = 1:pref.N; % Gives each firm an ID/name
color = linspecer(pref.N); % Generates a color for each of the firms that is evenly space out in / distinguishable.
sd = 0.5; % Standard deviation of each subpopulation
b = pref.boundary/2;


    %%% Population
    
    % Grid / Square
    [x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
    [X,Y] = meshgrid(x,y);

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
    weight = pref.n_ratio/(1+pref.n_ratio); % The left subpopulation share of total population
    F = F_l*weight + F_r*(1-weight); % Population pdf

    
    %%% Population descriptive statistics
    % Mean and covariance -- Johnson (1987)
    
    pop_mu = mu_l*weight + mu_r*(1-weight) % Population mean along each dimension; [mu_x mu_y]
    pop_sigma = sigma_l*weight + sigma_r*(1-weight) + weight*(1-weight)*(mu_l-mu_r)*(mu_l-mu_r)'; % Population covariance matrix
    pop_sd_x = sqrt(pop_sigma(1,1)) % Population std. dev. along x-axis; sigma_x


    %%% Firms
    % Randomy draw initial position of firms
    
    %xy_0 = (cos(rand(2,pref.N)*2*pi) .* (rand(2,pref.N) * 3))'; % start at (0,0) move uniformly up to 3 std. dev. in random direction
    xy_0 = normrnd(0,1,[2 pref.N])'; % bivariate normal distribution, mean 0, std. dev. 1, uncorrelated (rho=0)
    %xy_0 = (rand(2, pref.N)*pref.boundary-pref.boundary/2)'; % uniformly distributed on boundary

    
    %%% Preparing matrices for iteration
    
    % Creating a 3D matrix of firms. 
    % 1st dimension is firm, 2nd is (x,y) coordinates, and 3rd is iteration.
    xy = zeros(pref.N, 2, 2); 
    xy(:,:,1) = xy_0;
    
    % Creating a 3D matrix of market and utility.
    % 1st dimension is x, 2nd is y, and 3rd is iteration.
    market = zeros(pref.resolution, pref.resolution, 2); % The cell/value is the closest firm
    utility = zeros(pref.resolution, pref.resolution, 2); % The cell/value is utility
    %[market_0, utility_0, centroid_0, centroid_distance_0] = marketshare3(xy(:,:,1), [X(:) Y(:)]);
    [market_0, utility_0] = marketshare3(xy(:,:,1), [X(:) Y(:)]);
    market(:,:,1) = market_0;
    utility(:,:,1) = utility_0;
    
    % Number of customers for each firm
    F_firm = arrayfun(@(firm) ...
            sum(F(find(market_0==firm))), ...
            firms);
    
    % Market shares
    shares(1,:) = F_firm/sum(F(:)); % Calculate market share for each firm
    [~,rank(1,:)] = sort(shares(1,:), 'descend'); % Ranking firms according to market share
    
    % Eccentricity 
    % Firm euclidean distance from center/mean of voter distribution).
    eccentricity(1,:) = pdist2(xy(:,:,1), pop_mu, 'euclidean'); % in std. dev.
    mean_eccentricity(1,:) = mean(eccentricity(1,:));
    
    % Effective number of firms (ENP) -- Lever and Sergenti (2011), Laakso and Taagepera (1979)
    ENP(1,:) = sum(F(:))^2 / sum(F_firm(:).^2);
    % Herfindahl-Hirschman Index (HHI) -- Eiselt and Marianov (2011)
    % HHI(1,:) = sum( ( F_firm/sum(F(:)) ).^2 ); % 1/ENP
    
    % Misery
    % Quadratic loss function / representativeness
    misery(1,:) = sum(utility_0(:).*F(:))/sum(F(:));
  
    % Market perimeter
    market_props = regionpropsext(market(:,:,1), 'Extrema', 'Perimeter');
    perimeter(1,:) = cat(1,market_props.Perimeter)';
    perimeter_extrema(1,:) = cat(1,market_props.ExtremaPerimeter)';
    
    % Market centroid coordinates
    centroid = zeros(pref.N, 2, 2); 
    for firm=1:pref.N
        % Index for each customer of the firm.
        idx = find(market_0==firm);
        % Each xy-coordinat within the firm's market weighted with the probability density
        centroid(firm,:,1) = ([X(idx) Y(idx)]' * F(idx))' ./ F_firm(firm);
    end
    % Distance from each firm to the respective centroid of its market.
    centroid_distance(1,:) = diag(pdist2(xy(:,:,1), centroid(:,:,1), 'Euclidean'))';
    
%   centroid_distance(1,:)s = arrayfun(@(firm) ...
%             sum(centroid_distance_0(find(market_0==firm)) .* F(find(market_0==firm)))/F_firm(firm), ...
%             firms);
%   centroid_distance(1,:) = sum(centroid_distance_0(:).*F(:))/sum(F(:))
%   centroid_distance(1,:) = diag(pdist2(xy(:,:,1), cat(1,market_props.Centroid), 'Euclidean'))';


%% PLOTS

% 3D plot of density function
figure(10);
surf(x,y,F);%,'EdgeColor','none');
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-5 5 -5 5 0 .5])
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');

% 2D plot of contours
figure(11);
contour(x,y,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
hold on;
scatter( xy_0(:,1)' , xy_0(:,2)', [], color, 'filled'); % Plot the firms with respective colors.
xlim([-5 5]); ylim([-5 5]);
xlabel('x1'); ylabel('x2');
hold off;

bbox = [-b b b -b -b; b b -b -b b]'; % boundary box
[vx,vy] = voronoi(xy_0(:,1)' , xy_0(:,2)');
[V,C] = voronoin(xy_0);
% line intersect with boundary box
bbx = []; bby = []; %bbii = [];
for j = 2:length(vx)
    [xi, yi] = polyxpoly(bbox(:,1), bbox(:,2), vx(:,j), vy(:,j));
    bbx = [bbx; xi];
    bby = [bby; yi];
    %%bbii = [bbii; ii];
end
V3 = [bbx bby];

%patches = [-5 5; -5 1.2361; -2.0479 0.1954; -0.8703 0.9957; -2.5539 5.0000; -5 5];
%patches = [-0.8703 0.9957; 5.0000 1.7739; 5 -5; 0.8111 -5.0000; -2.2700 -1.6330; -2.0479 0.1954; -0.8703 0.9957];
%in = inpolygon(X,Y, patches(:,1)', patches(:,2)');

figure(12);
voronoi(xy_0(:,1)' , xy_0(:,2)');
xlim([-pref.boundary/2 pref.boundary/2]); ylim([-pref.boundary/2 pref.boundary/2]);
hold on;
%h = patch(patches(:,1)', patches(:,2)', 'red');
%set(h,'EdgeColor','none','FaceAlpha',0.25);
scatter(centroid(:,1,1)',centroid(:,2,1)', [], color, 'x');
hold off;

figure(14);
imagesc(utility_0)
set(gca,'ydir', 'normal');
cm=flipud(jet);
colormap(cm);

figure(13);
imagesc(market_0);
set(gca,'ydir', 'normal');
colormap(color); % Set the respective firm colors on each market share.
%hold on;
%scatter(centroid(:,1)',centroid(:,1)');
hold off;

