clearvars;


%% Setup

% Points - Random firm locations
n = 12;
%x0 = rand(1, n)*10-5; % Initial x-position of firm
%y0 = rand(1, n)*10-5; % Initial y-position of firm
%xy = [x0' y0'];
[x0, y0] = pol2cart( rand(n,1)*2*pi , rand(n,1)*3 );
xy = [x0 y0];

    
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


%% Plots

% Default plot
figure(1);
voronoi(xy(:,1)' , xy(:,2)', 'r');
%plot(vx,vy, 'b');
xlim([-6 6]); ylim([-6 6]);
%hold on;
%s = scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
%hold off;


DT = delaunayTriangulation(xy);
%tri = delaunay(xy(:,1)' , xy(:,2)');
figure(2);
triplot(DT, 'k')
xlim([-6 6]); ylim([-6 6]);

% figure(21);
% triplot(DT)
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% [centers, radii] = DT.circumcenter();
% theta = -pi:pi/20:pi;
% for iCircle=1:size(centers,1)
%   xCircle = centers(iCircle,1) + radii(iCircle)*cos(theta);
%   yCircle = centers(iCircle,2) + radii(iCircle)*sin(theta);
%   plot(xCircle, yCircle, 'k');
% end
% hold off;

xy_boundary = [xy; bbox(1:end-1,1) bbox(1:end-1,2)];
DT2 = delaunayTriangulation(xy_boundary);

triangles = size(DT2.ConnectivityList,1);

DT2X = NaN(size(DT2.ConnectivityList));
DT2Y = NaN(size(DT2.ConnectivityList));
for tri = 1:triangles
    DT2X(tri,:) = DT2.Points( DT2.ConnectivityList(tri,:), 1);
    DT2Y(tri,:) = DT2.Points( DT2.ConnectivityList(tri,:), 2);
end




%tri = delaunay(xy(:,1)' , xy(:,2)');
figure(3);
triplot(DT2, 'k')
xlim([-6 6]); ylim([-6 6]);
hold on;
for tri = 1:triangles
    patch(DT2X(tri,:), DT2Y(tri,:), tri, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
hold off;


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

centroid_tri_i(max_tri,:)


% l = sqrt(length(X(:)));
% DT_market_square = reshape(DT_market,[l l]);
% figure(5);
% imagesc(x,y,DT_market_square);
% set(gca,'ydir', 'normal');
% hold on;
% scatter( centroid_tri_i(:,1) , centroid_tri_i(:,2), 'r', 'filled');
% hold off;

% % Market centroid coordinates
% F_tri_i  = NaN(triangles,1);
% for tri= 1:triangles
%     % Index for each customer of the firm.
%     idx = find(market_i==firm);
% 
%     % Number of customers for each firm
%     F_tri_i(firm) = sum( F(idx) );
% 
%     % Each xy-coordinat within the firm's market weighted with the probability density
%     centroid(firm,:,i) = ([X(idx) Y(idx)]' * F(idx))' ./ F_tri_i(firm);
% end


%tri = delaunay(xy(:,1)' , xy(:,2)');

% figure(4);
% triplot(DT2)
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% [centers, radii] = DT2.circumcenter();
% theta = -pi:pi/20:pi;
% for iCircle=1:size(centers,1)
%   xCircle = centers(iCircle,1) + radii(iCircle)*cos(theta);
%   yCircle = centers(iCircle,2) + radii(iCircle)*sin(theta);
%   plot(xCircle, yCircle, 'k');
% end
% hold off;


xy_new = [xy; centroid_tri_i(max_tri,:)];

% Default plot
figure(7);
voronoi(xy_new(:,1)' , xy_new(:,2)', 'r');
%plot(vx,vy, 'b');
xlim([-6 6]); ylim([-6 6]);
%hold on;
%s = scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
%hold off;
