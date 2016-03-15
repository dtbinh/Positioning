% Population
pref.boundary = 10; % Number of standard deviations
pref.resolution = 50; % Length of the square (even number to include (0,0))
pref.mu = 1.5; % Mean of subpopulation
pref.n_ratio = 2; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation

sd = 0.5; % Standard deviation of each subpopulation
b = pref.boundary/2;
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
% The left subpopulation share of total population
weight = pref.n_ratio/(1+pref.n_ratio);
% Population probability density function (PDF). Left subpopulation is 
% left unchanged, while the right subpopulation will scale according to 
% n_ratio. Deviding by 2^2 has no effect on results, but scales PDF so 
% its comparable to the baseline distribution with mean (0,0) and 
% std. dev. (1,1). The subpopulations halves the std dev., thus their 
% variance will be one-fourth.
F = (F_l + F_r*1/pref.n_ratio)/4;
%F = F_l*weight + F_r*(1-weight); % Population pdf

% -----------

pref.repetitions = 200;
for rep=1:pref.repetitions
    
    pref.N = 6;
    [x0, y0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*3 );
    xy_test(1:pref.N,:,rep) = [x0 y0];
    %xy = [3.3735 0.7889; -0.1072 -3.4814; -3.9732 4.1955];
    xy_boundary = [xy_test(1:pref.N,:,rep); -5 5; 5 5; 5 -5; -5 -5];

    figure(20);
    clf reset; % Reset figure.
    [convexTriangleSets,DT] = triangles(xy_boundary);

    hulls = length(convexTriangleSets);
    F_tri  = NaN(hulls,1);
    centroid_tri  = NaN(hulls,2);
    share_tri  = NaN(hulls,1);
    for i=1:hulls
        points = DT.Points(DT(convexTriangleSets{i},:),:);
        k = convhull(points(:,1),points(:,2));

        % Index of customers within the convex hull
        idx = InPolygon(X(:),Y(:), points(k,1), points(k,2));

        % Number of customers in the convex hull
        F_tri(i) = sum( F(idx) );

        % The xy-coordinate of the centroid within the triangle weighted with the probability density
        centroid_tri(i,:) = ([X(idx) Y(idx)]' * F(idx))' ./ F_tri(i);


        [market_i, ~] = marketshare4([centroid_tri(i,:); xy_test(1:pref.N,:,rep)], [X(:) Y(:)]);
        idx2 = find(market_i==1);
        share_tri(i) = sum( F(idx2) ) / sum(F(:));


        %triangles = DT(convexTriangleSets{i},:)';
        %[~, idx] = unique(triangles);
        %triangles_s = triangles(sort(idx));
        %DT.Points(triangles_s,:)
    end

    [~, xy_new_idx] = max(F_tri);
    [share_max, share_max_idx] = max(share_tri);

    error_share(rep) = share_max-share_tri(xy_new_idx);
    error_distance(rep) = pdist([centroid_tri(xy_new_idx,:); centroid_tri(share_max_idx,:)], 'euclidean');

end

correct = (error_share == 0);
corrent_pct = sum(correct)./pref.repetitions;
[max_error_share, max_error_share_rep] = max(error_share);
max_error_distance = error_distance(max_error_share_rep);

mean_error_share = mean(error_share(~correct));
mean_error_distance = mean(error_distance(~correct));

summary = table(corrent_pct, max_error_share, mean_error_share, max_error_distance, mean_error_distance)



xy = xy_test(1:pref.N,:,max_error_share_rep);
figure(21);
DT = delaunayTriangulation([xy; -5 5; 5 5; 5 -5; -5 -5]);
figure(21);
triplot(DT, 'k');

