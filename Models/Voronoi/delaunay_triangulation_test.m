clearvars;


%% Setup

% Points - Random firm locations
n = 3;
x0 = rand(1, n)*10-5; % Initial x-position of firm
y0 = rand(1, n)*10-5; % Initial y-position of firm
xy = [x0' y0'];
%[x0, y0] = pol2cart( rand(n,1)*2*pi , rand(n,1)*3 );
%xy = [x0 y0];

    
%% Boundary box

% Boundary box
%bbox = [x([1 end end 1 1]); y([end end 1 1 end])]';
bbox = [-5 5 5 -5 -5; 5 5 -5 -5 5]';

% Boundary box line segments
bvx = [ bbox(1:end-1,1)'; circshift( bbox(1:end-1,1)', [0,-1] ) ];
bvy = [ bbox(1:end-1,2)'; circshift( bbox(1:end-1,2)', [0,-1] ) ];


%% Population
pref.boundary = 10; % Number of standard deviations
pref.resolution = 50; % Length of the square (even number to include (0,0))
pref.mu = 0; % Mean of subpopulation
pref.n_ratio = 1; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation

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


%% Calculations

% Delaunay Triangulation without boundaries.
DT = delaunayTriangulation(xy);
% Delaunay Triangulation with boundaries.
xy_boundary = [xy; bbox(1:end-1,1) bbox(1:end-1,2)];
DT2 = delaunayTriangulation(xy_boundary);

% Number of triangles
triangles = size(DT2.ConnectivityList,1);

% Get coordinates for each triangle
DT2X = NaN(size(DT2.ConnectivityList));
DT2Y = NaN(size(DT2.ConnectivityList));

%perimeter = NaN(triangles,1);
%area = NaN(triangles,1);
for tri = 1:triangles
    DT2points = DT2.Points( DT2.ConnectivityList(tri,:), :);
    DT2X(tri,:) = DT2points(:, 1);
    DT2Y(tri,:) = DT2points(:, 2);
    
    %length_sides = sqrt( sum( diff([DT2points; DT2points(1,:)], [], 1).^2, 2) );
    %perimeter(tri) = sum(length_sides);
    
    % Heron's Formula
    %s = perimeter(tri)/2;
    %area(tri) = sqrt(s * prod(s-length_sides));
end

% Calculate the length of each side in triangle
length_sides = sqrt( diff([DT2X DT2X(:,1)], [], 2).^2 + diff([DT2Y DT2Y(:,1)], [], 2).^2 );
perimeter = sum(length_sides, 2);
% Semi-perimeter
s = perimeter/2;
% Heron's Formula for area
area = sqrt(s .* prod(repmat(s,1,3)-length_sides, 2));
% Radius of the incirlce
inradius = area ./ s;

% Alternative eccentricity measure
eccentricity = perimeter./sqrt(area);


% Calculate the market share of each triangle and the centroid
DT_market = NaN( length(X(:)), 1);
F_tri_i  = NaN(triangles,1);
centroid_tri_i  = NaN(triangles,2);
for tri = 1:triangles
    idx = inpolygon(X(:), Y(:), DT2X(tri,:), DT2Y(tri,:));
    DT_market(idx) = tri;
    F_tri_i(tri) = sum( F(idx) );
    % Each xy-coordinat within the triangle's market weighted with the probability density
    centroid_tri_i(tri,:) = ([X(idx) Y(idx)]' * F(idx))' ./ F_tri_i(tri);
end

% Market share of each triangle.
%F_tri_i/sum(F(:))
[~, max_tri] = max(F_tri_i/sum(F(:)));
[~, max_tri2] = max( F_tri_i/sum(F(:)) .* eccentricity );
[~, max_tri3] = max( F_tri_i/sum(F(:)) .* log(1+1./inradius) );

% The optimal location
centroid_tri_i(max_tri,:);

% Color triangle based on customer share.
alpha_tri = F_tri_i./sum(F_tri_i);

% firm coordinates including new firm
xy_new = [xy; centroid_tri_i(max_tri,:)];


%% Plots

% Voronoi plot
figure(1);
voronoi(xy(:,1)' , xy(:,2)', 'k');
xlim([-5 5]); ylim([-5 5]);
set(gca,'YTick',(-5:5:5));
set(gca,'XTick',(-5:5:5));
hold on;
s = scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
hold off;

% Delaunay Triangulation with boundaries
figure(3);
triplot(DT2, ':k')
xlim([-5 5]); ylim([-5 5]);
set(gca,'YTick',(-5:5:5));
set(gca,'XTick',(-5:5:5));
hold on;
triplot(DT, '--k'); % 'color', [0.5 0.5 0.5]
scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
scatter( centroid_tri_i(max_tri,1)' , centroid_tri_i(max_tri,2)', 'xr');
scatter( centroid_tri_i(max_tri2,1)' , centroid_tri_i(max_tri2,2)', '+g');
scatter( centroid_tri_i(max_tri3,1)' , centroid_tri_i(max_tri3,2)', '+c');
%for tri = 1:triangles
    %patch(DT2X(tri,:), DT2Y(tri,:), tri, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', alpha_tri(tri), 'EdgeColor', 'none');
%    patch(DT2X(tri,:), DT2Y(tri,:), tri, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', eccentricity(tri)/max(eccentricity), 'EdgeColor', 'none');
%end
hold off;


% Voronoi plot Wtih new firm
figure(7);
voronoi(xy_new(:,1)' , xy_new(:,2)', 'k');
xlim([-5 5]); ylim([-5 5]);
set(gca,'YTick',(-5:5:5));
set(gca,'XTick',(-5:5:5));
hold on;
scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
scatter( centroid_tri_i(max_tri,1)' , centroid_tri_i(max_tri,2)', 'filled', 'r');
hold off;
