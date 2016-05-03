% TO-DO: Create population function with input [x y pref] and output F or [F, F_l, F_r]


%% Population distribution
pref.boundary = 10; % Number of standard deviations
pref.resolution = 200; % Length of the square (even number to include (0,0))
b = pref.boundary/2;
[x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
[X,Y] = meshgrid(x,y);
%Population
pref.mu = 1.5; % Mean of subpopulation
pref.n_ratio = 2; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
sd = 0.5; % Standard deviation of each subpopulation
mu_r = [pref.mu 0]; % Subpopulation mean only deviate on x-axis
sigma_r = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r); % Subpopulation pdf evaluated at each point in grid/square
F_r = reshape(F_r,length(y),length(x)); % Formated as grid
mu_l = [-pref.mu 0]; % Subpopulation mean only deviate on x-axis
sigma_l = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l); % Subpopulation pdf evaluated at each point in grid/square
F_l = reshape(F_l,length(y),length(x)); % Formated as grid
weight = pref.n_ratio/(1+pref.n_ratio);
F = (F_l + F_r*1/pref.n_ratio)/4;


%% Loation of firms

pref.N = 7;
[x0, y0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*5 );
xy = [x0 y0];
color = brewermap(pref.N,'Paired');
%xy1=xy
%xy = xy1
%pref.N = 4;
%color = brewermap(pref.N,'Paired');

%% Procedure

% CONVEX HULL OF ALL POINTS
convh = convhull(xy(:,1)' , xy(:,2)');

% VORONOI AND DELAUNAY
[vx,vy] = voronoi(xy(:,1)' , xy(:,2)');
vexy = [vx(1,:)' vy(1,:)' vx(2,:)' vy(2,:)']; % V edge as [x1 y1 x2 y2]

DT = delaunayTriangulation(xy(:,1), xy(:,2));
dexy = reshape(DT.Points(DT.edges',:)', 4, [])'; % D edge as [x1 y1 x2 y2]

% VORONOI VERTICES (WHERE THREE OR MORE EDGES MEET)
[vxy_all, ~, vxy_i] = unique( reshape([vx; vy]', [], 2), 'rows'); % Returns sorted array
n = histc(vxy_i, unique(vxy_i));
vxy = vxy_all(vxy_i(find(n >= 3)),:);
% Density mass at Voronoi verticies
F_v3xy = (mvnpdf(vxy, mu_l, sigma_l) + mvnpdf(vxy, mu_r, sigma_r)*1/pref.n_ratio)/4;


% INTERSECTION POINTS BETWEEN VORONOI EDGES AND DELAUNAY EDGES
% Each generating points is a maximum point in the cone landscape, the 
% Voronoi vertices (where 3 or more edges meet) are (local) minimum points. 
% While the intersection points are saddle points in the cone landscape. 
intersectvd = lineSegmentIntersect( vexy, dexy );
% Ignore Delaunay edges with multiple intersection points, since these
% intersection points are not on the surface of the cone landscape.
intersectdsingle = (sum(intersectvd.intAdjacencyMatrix, 1)<2);
% Following variables use the single Delaunay edge intersection condition.
% Adjeceny matrix of Voronoi and Delaunay edges intersection
intersectvdadj = logical(intersectvd.intAdjacencyMatrix .* repmat(intersectdsingle, size(vexy,1), 1) );
% Coordinates of intersection points
intersectvdxy = [intersectvd.intMatrixX(intersectvdadj) intersectvd.intMatrixY(intersectvdadj)];
% Index of the voronoi edges (with single Delaunay edge intersection condition)
%intersectvsingle = find(sum(intersectvdadj, 2));
% Density mass at intersection points
F_intersectvdxy = (mvnpdf(intersectvdxy, mu_l, sigma_l) + mvnpdf(intersectvdxy, mu_r, sigma_r)*1/pref.n_ratio)/4;

% HALF DISTANCES
% Calculate the half distance between pairs of points and the combination 
% of pairs (Only considering pairs where direct line between points cross a 
% single Voronoi line, ie. single Delaunay edge intersection condition).
[xys_hdist, combis] = pdist2(xy, intersectvdxy, 'euclidean', 'Smallest', 2);
xys_hdist(2,:) = []; % Remove redundant distances (equal distance to Voronoi edges) 
combis = sort(combis,1)'; % Transepose and sort combinations.
[xys_hdist_sort, sort_i] = sort(xys_hdist, 'descend'); % Sorted with largest distance first

F_intersectvdxy_sort = F_intersectvdxy(sort_i);


% % Index for all combination of pair of points
% combi = nchoosek(1:pref.N, 2); % fliplr() flip so order matches pdist
% % Distance between each pair of points sorted with largest distance first
% xy_dist = pdist(xy, 'euclidean');
% %xy_dist = diag(pdist2(xy(combi(:,1),:), xy(combi(:,2),:), 'euclidean'));
% xy_dist_sort = sort(xy_dist, 'descend');


% Margin of error
err1 = 0.0000000001;
err2 = 1.05;

% LOOP THROUGH ALL HALF DISTANCCES SETTING RADIUS
wells = cell(length(xy_dist_sort),1);
%wellx = cell(length(xy_dist_sort),1);
%welly = cell(length(xy_dist_sort),1);
intersectxy_inhullx = cell(length(xy_dist_sort),1);
intersectxy_inhully = cell(length(xy_dist_sort),1);

for vd = sort_i
    % Radius is half distance between pair of points.
    %r = xy_dist_sort(5)/2; %1.4; % 1.7
    %r = xys_hdist_sort(vd); %1.4; % 1.7
    r = xys_hdist(vd); %1.4; % 1.7

    
    % Find Voronoi verticies outside circles / within well.
    vxy_dist = pdist2(xy, vxy, 'euclidean', 'Smallest', 1);
    vxy_active = find((vxy_dist > r));

    
    % Intersection points between the circles / first type of polygon points for the well
    % For a given r (radius of cirlce) strictly consider the cirlces with two intersection points.
    %combi_active = find((xy_dist < r*2));
    combi_active = find((xys_hdist <= r));
    intersectcirclexy = nan(1,2);
    for c=combi_active
        %c1 = combi(c,1);
        %c2 = combi(c,2);
        c1 = combis(c,1);
        c2 = combis(c,2);
        %[xy(c1,1) xy(c1,2) r xy(c2,1) xy(c2,2) r]
        [xout, yout] = circcirc(xy(c1,1), xy(c1,2), r+err1, xy(c2,1), xy(c2,2), r);
        intersectcirclexy(end+1:end+2,:) = [xout' yout'];
    end
    %intersectcirclexy = intersectcirclexy(~any(isnan(intersectcirclexy),2),:); % Delete all [nan nan] rows (first and any missing intersection).
    intersectcirclexy(1,:) = []; % Delete first [nan nan] row.
    %intersectcirclexy = unique(intersectcirclexy, 'Rows'); % Delete duplicates

    
    % Find circle intersection points with the convex hull of all points / second type of polygon points for the well
    conhxy = [xy(convh(1:end-1),:) circshift(xy(convh(1:end-1),:),-1)]; % Convex hull as [x1 y1 x2 y2].
    intersecthullxy_all = linecirc2( conhxy, [xy repmat(r, pref.N, 1)] );
    %convh_vector = diff(xy(convh,:)); % A vector for each face of the polygon
    %convh_norm = sqrt(dot(convh_vector, convh_vector, 2)); % Norm of vector
    %convh_normvector = convh_vector ./ repmat(convh_norm, 1, 2); % Normalise vector to length 1.
    %% Intersection point is the vector with length r (radius of circle) going both CW and CCW around the polygon. 
    %intersecthullxy_all = [xy(convh(1:end-1), :) + convh_normvector * r; ...
    %                   xy(flipud(convh(1:end-1)), :) - circshift(flipud(convh_normvector),-1) * r];
    %% Remove the convex hull intersection points that fall within one of the circles. 
    intersecthull_seperate = (pdist2(xy, intersecthullxy_all, 'euclidean', 'Smallest', 1) >= r-err1)';
    intersecthullxy = intersecthullxy_all(intersecthull_seperate,:);
    
    % Ignore the intersection points outside of the convex hull of all points
    intersectxy_all = [intersectcirclexy; intersecthullxy];
    inhull = InPolygon(intersectxy_all(:,1), intersectxy_all(:,2), xy(convh,1) , xy(convh,2));
    intersectxy_inhull = intersectxy_all(inhull,:);

    
    % TO-DO: If circles fully conver the convex hull of points. Ie. no wells 
    % TO-DO: Multiple wells -- Sorts points CW around each V vertex.
    ii = (pdist2(xy, intersectvdxy, 'euclidean', 'Smallest', 1) <= r);
    %intersectvdxy(ii, :)
    % D lines with single intersection points
    dexy_single = dexy(intersectdsingle',:);
    
    % Calculate all convex hull around points that do no intersect D lines
    % with single intersection point.
    ch = convhullc(intersectxy_inhull, dexy_single(ii,:));
    wells{vd} = ch;
    % D lines with single intersection point and potential seperating midpoint
    %dexy_single_temp = reshape([dexy_single(ii',:) nan(sum(ii),2)]', 2, [])'
    %
    %well_r = convhull( intersectxy_inhull(:,1), intersectxy_inhull(:,2) );

    %wellx{vd} = intersectxy_inhull(well_r,1);
    %welly{vd} = intersectxy_inhull(well_r,2);
    intersectxy_inhullx{vd} = intersectxy_inhull(:,1);
    intersectxy_inhully{vd} = intersectxy_inhull(:,2);
end




% fig10 = figure(10);
% clf(fig10);
% %imagesc(x, y, F);
% contour(x,y,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
% hold on;
% for n=1:pref.N
%     rectangle('Position',[xy(n,1)-r xy(n,2)-r 2*r 2*r],'Curvature', 1, 'EdgeColor', 'none', 'FaceColor', [0.9 0.9 0.9])
%     %plot(xy(n,1)+xp, xy(n,2)+yp, 'Color', color(n,:));
%     %plot(xy(n,1), xy(n,2), '.r', 'MarkerSize', r)
% end
% scatter( xy(:,1) , xy(:,2), [], color, 'filled');
% text(xy(:,1), xy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','black');
% plot( xy(convh,1) , xy(convh,2),'m-');
% scatter( intersectxy_all(inhull,1), intersectxy_all(inhull,2), [], 'b');
% plot( intersectxy_inhull(well_r,1) , intersectxy_inhull(well_r,2), 'Color', [255 165 0]/255);
% %scatter( intersectcirclexy(:,1), intersectcirclexy(:,2), [], 'r');
% scatter( intersecthullxy_all(:,1), intersecthullxy_all(:,2), [], 'k', 'filled');
% %scatter( intersecthullxy2(:,1), intersecthullxy2(:,2), [], 'b');
% %text(intersecthullxy(:,1), intersecthullxy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','blue');
% xlim([-5 5]); ylim([-5 5]);
% hold off;
% 
% 
% figure(11);
% scatter( xy(:,1) , xy(:,2), [], 'k', 'filled');
% xlim([-5 5]); ylim([-5 5]);
% hold on;
% plot(vx,vy, 'b');
% triplot(DT, ':', 'Color', [0.5 0.5 0.5]);
% scatter( intersectvdxy(:,1) , intersectvdxy(:,2), [], 'r');
% scatter( 0, 0, [], 'b*');
% plot( intersectxy_inhull(well_r,1) , intersectxy_inhull(well_r,2), 'Color', [255 165 0]/255);
% scatter( vxy(vxy_active,1), vxy(vxy_active,2), 'r', 'filled'); %F_v3xy(vxy_active)/max(F(:))*100*36
% hold off;


for i=sort_i
%i=3;

figure(11);
scatter( xy(:,1) , xy(:,2), [], 'k', 'filled');
xlim([-5 5]); ylim([-5 5]);
hold on;
plot(vx,vy, 'b');
triplot(DT, ':', 'Color', [0.5 0.5 0.5]);
scatter( intersectxy_inhullx{i} , intersectxy_inhully{i}, [], 'r');
scatter( 0, 0, [], 'b*');
%plot( wellx{i} , welly{i}, 'Color', [255 165 0]/255);
for w=1:length(wells{i})
   plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [255 165 0]/255); 
end
hold off;
colorbar;

fig10 = figure(10);
r =  xys_hdist(i);
clf(fig10);
contour(x,y,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
hold on;
for n=1:pref.N
    rectangle('Position',[xy(n,1)-r xy(n,2)-r 2*r 2*r],'Curvature', 1, 'EdgeColor', 'none', 'FaceColor', [0.9 0.9 0.9])
end
scatter( xy(:,1) , xy(:,2), [], color, 'filled');
text(xy(:,1), xy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','black');
plot( xy(convh,1) , xy(convh,2),'m-');
%plot( wellx{i} , welly{i}, 'Color', [255 165 0]/255);
%scatter( intersecthullxy(:,1), intersecthullxy(:,2), 'g');
scatter( intersectxy_inhullx{i} , intersectxy_inhully{i}, [], 'k', 'filled');
%scatter(intersectvdxy(ii,1), intersectvdxy(ii,2), 20, 'b', 'filled');
%plot( dexy_single_temp(:,1), dexy_single_temp(:,2) , 'b');
xlim([-5 5]); ylim([-5 5]);
for w=1:length(wells{i})
   plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [255 165 0]/255); 
end
hold off;
colorbar;

pause(0.5)
end




% u(1,1)
% F(1,1)
% r = xys_hdist_sort(i);
% d = F_intersectvdxy_sort(i);
% %[ u(25,25) F(25,25) r d  abs(u(25,25)-r) F(25,25)-d  abs(u(25,25)-r)/(F(25,25)-d)]
% [ u(25,25) F(25,25) ...
%     r d ...
%     u(25,25)/r F(25,25)/d ...
%     (F(25,25)/d > u(25,25)/r) ...
%     abs(u(25,25)-r) / max(F(:))/F(25,25) ...
%     r / max(F(:))/F(25,25) ...
%     F(25,25) * r/u(25,25) ... 
%     F(25,25)/u(25,25) * r]
% 
% min(u(:))


color = brewermap(pref.N,'Paired');
%[~, ii] = sort(xys_hdist ./ F_intersectvdxy', 'descend');
[F_max, maxi] = max(F(:));
xy_max = [X(maxi) Y(maxi)];
F2 = F / F_max;
u = reshape(pdist2(xy, [X(:) Y(:)], 'euclidean', 'Smallest', 1), size(X));
u2 = u / max(u(:));
[~, sort_i2] = sort( (F_intersectvdxy/F_max)' ./ (xys_hdist / max(u(:))), 'descend' );

for i=sort_i2
%i=2;
r =  xys_hdist(i);
r2 = r / max(u(:));
r3 = pdist([intersectvdxy(i,:); xy_max], 'euclidean');
F_i = F_intersectvdxy(i);
F2_i = F_i / F_max;
%u4 = (F ./ u >  F_intersectvdxy(i) / r);
%Fabove = (F2 > F_intersectvdxy(i)/F_max);
Fabove = (F2/F2_i > u2/r2);
%F_new(Fabove) = F2(Fabove);
%distbelow = (u(Fabove) < r);
%F_new(~distbelow) = F2(~distbelow);
%F_new(distbelow) = F2(distbelow) - F2(distbelow) .* u(distbelow)/max(u(distbelow));
%F_new(Fabove) = F2(Fabove) .* ( 1 + (u(Fabove)-r) ./ u(Fabove) ) ;
%u5(u4) = F(u4) - (F(u4) - F_intersectvdxy(i) * (u(u4) - r) );  % .* ( 1 + (u(u4)-r) ./ u(u4) ) ;
%F_new(~Fabove) = 0;
%F_new = reshape(F_new, size(F2));
%F_new = F2 - (1-F2_i)/r3;
F_new = F2 - (F2 - u2) * (F2_i/r2) ;
%F_new(Fabove) = F2(Fabove) - (F2(Fabove) - u2(Fabove)) * (F2_i/r2) ;
%F_new(~Fabove) = F2(~Fabove) + (F2(~Fabove) - u2(~Fabove)) * (F2_i/r2) ;
fig12 = figure(12);
clf(fig12);
hold on;
for n=1:pref.N
    rectangle('Position',[xy(n,1)-r xy(n,2)-r 2*r 2*r],'Curvature', 1, 'EdgeColor', 'none', 'FaceColor', [0.9 0.9 0.9])
end
h = imagesc(x, y, F_new )   ;
set(h, 'AlphaData', Fabove);
colorbar;
caxis([-1,1]);
colormap(fig12,'default')
map = colormap;
map2(33:64,:) = map(1:2:end,:);
map2(32,:) = [0 0 0];
map2(1:31,:) = repmat([0.7 0.7 0.7], 31, 1);
colormap(map2)
%contour(x,y,F,[F_intersectvdxy(2):0.03:max(F(:))]);
%plot( xy(convh,1) , xy(convh,2),'m-');
scatter( intersectvdxy(i,1) , intersectvdxy(i,2), 'filled', 'MarkerFaceColor', 'm');
for w=1:length(wells{i})
   plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [255 165 0]/255); 
end
%scatter( xy(:,1) , xy(:,2), 50, color, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 2);
contour(x, y, u2, [xys_hdist(i)/max(u(:)) xys_hdist(i)/max(u(:))], 'Color', 'black');
xlim([-5 5]); ylim([-5 5]);
set(gca,'ydir', 'normal');
hold off;
pause(1)
end




figure(13);
imagesc(x, y, u);
hold on;
scatter( intersectvdxy(i,1) , intersectvdxy(i,2), 'k*');
contour(x, y, u2, [xys_hdist(i)/max(u(:)) xys_hdist(i)/max(u(:))], 'Color', 'black');
hold off;
xlim([-5 5]); ylim([-5 5]);
set(gca,'ydir', 'normal');
colormap(flipud(winter(25)));
colorbar;


figure(14);
%imagesc(x, y, F);
imagesc(x, y, F/max(F(:)));
colorbar;
xlim([-5 5]); ylim([-5 5]);
set(gca,'ydir', 'normal');




%% TEMPORARY


% Vectors to draw circle
%ang=0:0.01:2*pi; 
%xp=r*cos(ang);
%yp=r*sin(ang);



% Finding extrema points of distribution
extremax = extrema(pref.mu, pref.n_ratio); % [x_saddle x_right x_left]
extremaxy = [extremax' extremax'*0]; % Extrema y coordinate is zero.
% Finding density mass at extrema
extremaF = (mvnpdf(extremaxy, mu_l, sigma_l) + mvnpdf(extremaxy, mu_r, sigma_r)*1/pref.n_ratio)/4;




%DT.Points
%unique( reshape(dexy(intersectdsingle,:)', 2, [])', 'rows');

% Generating points (with single Delaunay edge intersection condition)
%xysingle = unique( reshape([vx(:,intersectvsingle); vy(:,intersectvsingle)]', [], 2), 'rows');

%intersect(xy, xysingle, 'rows')


%[~, combi_order] = sort(combi(:,1));
%combi = combi(combi_order,:)




unique(vxy_i)

%vxy_i(2) = []
%vxy_i = [vxy_i(1:6); 3; vxy_i(7:end)]
% Find index of first coordinate that have at least three identical arrays
%([diff(vxy_i,2);0.1]==0) .* (diff(vxy_i,1)==0)

find( (diff(vxy_i,1)==0).*([1; diff(diff(vxy_i,1))] ~= 0) )
vxy_i3 = circshift(vxy_i, -2);
vxy_i3(end-1:end) = 0;
vxy_i-vxy_i3
[ (diff(vxy_i,1)==0) ([1; diff(diff(vxy_i,1))] ~= 0) ]
%vxy_i = find(diff(vxy_i,2) == 0);



%imagesc(x, y, u);
%imagesc(x, y, (u > r) );
%imagesc(x, y, F );
%imagesc(x, y, (F > F_intersectvdxy(1)) );
%imagesc(x, y, (F > u2) );



%text(xy(:,1), xy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','black');
%scatter( intersectxy_all(inhull,1), intersectxy_all(inhull,2), [], 'b');
%plot( intersectxy_inhull(well,1) , intersectxy_inhull(well,2), 'Color', [255 165 0]/255);
%scatter( intersectcirclexy(:,1), intersectcirclexy(:,2), [], 'r');
%scatter( intersecthullxy(:,1), intersecthullxy(:,2), [], 'k', 'filled');
%scatter( intersecthullxy2(:,1), intersecthullxy2(:,2), [], 'b');
%text(intersecthullxy(:,1), intersecthullxy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','blue');


%u2 = u * ( max(F(:)) - F_intersectvdxy(2) );
%u3 = (F ./ ( max(F(:)) - F_intersectvdxy(2) ) > u) ; % (F > F_intersectvdxy(3));% .* 
%u5(u4) = F(u4) - F(u4) ./ u(u4) * xys_hdist_sort(i) ;
% .* F; %(F / F_intersectvdxy_sort(2) >  u / xys_hdist_sort(2)) .* F;