% Population
% version 0.01
% Jonas K. Sekamane

% Inspired in part by: 
%   Lever and Sergenti (2011)

clearvars;

rng('default'); % Seed such that the randomly generated results are repeatable

% Preferences
pref.boundary = 10; % Number of standard deviations
pref.resolution = 70; % Length of the square (even number to include (0,0))
pref.N = 5; % Number of firms
pref.mu = 1.5; % mean of subpopulation
pref.n_ratio = 2; % relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation

% Setup
unit = pref.resolution/pref.boundary;rng('default'); % Seed such that the randomly generated results are repeatable

% Preferences
pref.boundary = 10; % Number of standard deviations
pref.resolution = 70; % Length of the square (even number to include (0,0))
pref.N = 5; % Number of firms
pref.mu = 1.5; % mean of subpopulation
pref.n_ratio = 2; % relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation

% Setup
unit = pref.resolution/pref.boundary; % 1 standard deviation
firms = 1:pref.N; % Gives each firm an ID/name
color = linspecer(pref.N); % Generates a color for each of the firms that is evenly space out in / distinguishable.
sd = 0.5; % standard deviation of each subpopulation
b = pref.boundary/2;

% Grid / Square
[x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
[X,Y] = meshgrid(x,y);

% Subpopulation right
mu_r = [pref.mu 0]; % subpopulation mean only deviate on x-axis
sigma_r = [sd^2 0; 0 sd^2]; % subpopulation sd. Uncorrelated bivariate distribution, rho=0;
F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r);
F_r = reshape(F_r,length(y),length(x));

% Subpopulation left
mu_l = [-pref.mu 0];
sigma_l = [sd^2 0; 0 sd^2];
F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l);
F_l = reshape(F_l,length(y),length(x));

% Total population 
weight = pref.n_ratio/(1+pref.n_ratio);
F = F_l*weight + F_r*(1-weight);
% Population mean and covariance ? Johnson (1987)
pop_mu = mu_l*weight + mu_r*(1-weight)
pop_sigma = sigma_l*weight + sigma_r*(1-weight) + weight*(1-weight)*(mu_l-mu_r)*(mu_l-mu_r)';
pop_sd_x = sqrt(pop_sigma(1,1))


% Randomy draw initial position of firms
%xy0 = (cos(rand(2,pref.N)*2*pi) .* (rand(2,pref.N) * 3))'; % start at (0,0) move uniformly up to 3 std. dev. in random direction
xy0 = normrnd(0,1,[2 pref.N])'; % bivariate normal distribution, mean 0, std. dev. 1, uncorrelated (rho=0)
%xy0 = (rand(2, pref.N)*pref.boundary-pref.boundary/2)'; % uniformly distributed on boundary

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
scatter( xy0(:,1)' , xy0(:,2)', [], color, 'filled'); % Plot the firms with respective colors.
xlim([-5 5]); ylim([-5 5]);
xlabel('x1'); ylabel('x2');
hold off;

bbox = [-b b b -b -b; b b -b -b b]'; % boundary box

[vx,vy] = voronoi(xy0(:,1)' , xy0(:,2)');
[V,C] = voronoin(xy0);
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
voronoi(xy0(:,1)' , xy0(:,2)');
xlim([-pref.boundary/2 pref.boundary/2]); ylim([-pref.boundary/2 pref.boundary/2]);
hold on;
%h = patch(patches(:,1)', patches(:,2)', 'red');
%set(h,'EdgeColor','none','FaceAlpha',0.25);
hold off;

[market, utility] = marketshare3(xy0, [X(:) Y(:)]);
F_firm = arrayfun(@(firm) sum(F(find(market==firm))), firms);
shares = F_firm/sum(F(:)); % Calculate share for each firm

[~,rank] = sort(shares, 'descend');

% Eccentricity (euclidean from center/mean of voter distribution).
eccentricity = pdist2(xy0, pop_mu, 'euclidean'); % in std. dev.
%eccentricity = sqrt( sum( (xy0-repmat(pop_mu, pref.N,1)).^2 ,2))/10; 
mean(eccentricity);
%mean(pdist2([X(:) Y(:)],pop_mu,'euclidean')); % eccentricity population

% Effective number of firms (ENP) ? Lever and Sergenti (2011), Laakso and Taagepera (1979)
ENP = sum(F(:))^2 / sum(F_firm(:).^2);
% Herfindahl-Hirschman Index (HHI) ? Eiselt and Marianov (2011)
% HHI = sum( ( F_firm/sum(F(:)) ).^2 ); % 1/ENP

% Misery (quadratic loss function)
%misery = arrayfun(@(firm) sum( pdist2([X(find(market==firm)) Y(find(market==firm))], xy0(firm,:), 'euclidean') .* F(find(market==firm)) )/100, firms);
%misery = arrayfun(@(firm) sum( utility(find(market==firm)).* F(find(market==firm)) )/100, firms);
misery = sum(utility(:).*F(:))/sum(F(:));

figure(14);
imagesc(utility)
set(gca,'ydir', 'normal');
cm=flipud(jet);
colormap(cm);

figure(13);
imagesc(market);
set(gca,'ydir', 'normal');
colormap(color); % Set the respective firm colors on each market share.
hold on;


