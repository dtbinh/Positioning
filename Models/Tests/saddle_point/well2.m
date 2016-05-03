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

pref.N = 12;
[x0, y0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*5 );
xy = [x0 y0];
color = brewermap(pref.N,'Paired');


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
% Delaunay edges with single intersection points
dexysingle = dexy(intersectdsingle',:);
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
[xys_hdist_sort, sort_dist] = sort(xys_hdist, 'descend'); % Sorted with largest distance first


% Margin of error
err1 = 0.0001; % 0.0000000001;


% LOOP THROUGH ALL HALF DISTANCCES SETTING RADIUS
wells = cell(length(xys_hdist),1);
points_inhull = cell(length(xys_hdist),1);

for vd = sort_dist
    % Radius is half distance between pair of points.
    r = xys_hdist(vd);

    % Find Voronoi verticies outside circles / within well.
    vxy_dist = pdist2(xy, vxy, 'euclidean', 'Smallest', 1);
    vxy_active = find((vxy_dist > r));
    
    % FIRST TYPE OF POTENTIAL POLYGON POINTS FOR THE WELLS
    % Intersection points between the circles
    % For a given r (radius of cirlce) strictly consider the cirlces with two intersection points.
    combi_active = find((xys_hdist <= r));
    intersectcirclexy = nan(1,2);
    for c=combi_active
        c1 = combis(c,1);
        c2 = combis(c,2);
        [xout, yout] = circcirc(xy(c1,1), xy(c1,2), r+err1, xy(c2,1), xy(c2,2), r);
        intersectcirclexy(end+1:end+2,:) = [xout' yout'];
    end
    %intersectcirclexy = intersectcirclexy(~any(isnan(intersectcirclexy),2),:); % Delete all [nan nan] rows (first and any missing intersection).
    intersectcirclexy(1,:) = []; % Delete first [nan nan] row.
    %intersectcirclexy = unique(intersectcirclexy, 'Rows'); % Delete duplicates

    % SECOND TYPE OF POTENTIAL POLYGON POINTS FOR THE WELLS
    % Find circle intersection points with the convex hull of all points
    conhxy = [xy(convh(1:end-1),:) circshift(xy(convh(1:end-1),:),-1)]; % Convex hull as [x1 y1 x2 y2].
    intersecthullxy_all = linecirc2( conhxy, [xy repmat(r, pref.N, 1)] );
    
    % COMBINING THE TWO TYPES OF POTENTIAL POLYGON POINTS
    % Ignore the intersection points outside of the convex hull of all points
    intersectxy_all = [intersectcirclexy; intersecthullxy_all];
    % Remove the intersection points that fall within one of the circles. 
    intersectxy_seperate = (pdist2(xy, intersectxy_all, 'euclidean', 'Smallest', 1) >= r-err1)'; %>= r-err1
    intersect_all = intersectxy_all(intersectxy_seperate,:);
    inhull = InPolygon(intersect_all(:,1), intersect_all(:,2), xy(convh,1) , xy(convh,2));
    intersect_inhull = intersect_all(inhull,:);

    % TO-DO: If circles fully conver the convex hull of points. Ie. no wells 

    % SPLITTING POTENTIAL POLYGON POINTS INTO GROUPS & CREATING CONVEX HULL 
    % Identify the intersection points between Delaunay edges and Voronoi 
    % edges that fall within the radius of the circle. The intersection 
    % point or the Delaunay edges on which they are place divide the points 
    % into different wells.
    intersectvd_incircle = (pdist2(xy, intersectvdxy, 'euclidean', 'Smallest', 1) <= r);
    % Find all wells/convex hull around points that do no intersect 
    % specific Delaunay edges (single Delaunay edge intersection condition,
    % and intersection point within circle radius).
    ch = convhullc(intersect_inhull, dexysingle(intersectvd_incircle,:));
    
    % TO-DO: disregard wells that cover more than half the circle/where 
    % generating point falls within the well.
    for chj=1:length(ch)
        inwell = InPolygon(xy(:,1), xy(:,2), ch{chj}(:,1) , ch{chj}(:,2));
        if any(inwell)
            ch{chj} = [];
        end
    end
    ch = ch(~cellfun('isempty',ch));
    
    % Store the coordinates the every well
    wells{vd} = ch;
    % Store all potential polygon points
    points_inhull{vd} = intersect_inhull;
    
end



%% PLOTS


for i=sort_dist
    %TO-DO: Skip first, since it will not create a well, but a point/line.
    %TO-DO: Stop once no well exists.
    %i=2;
    r =  xys_hdist(i);

    figure(11);
    scatter( xy(:,1) , xy(:,2), [], 'b', 'filled');
    hold on;
    plot(vx,vy, 'b');
    triplot( DT, ':', 'Color', [0.5 0.5 0.5]);
    scatter( points_inhull{i}(:,1) , points_inhull{i}(:,2), [], 'k', 'filled');
    scatter( 0, 0, [], 'r*');
    for w=1:length(wells{i})
       plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [255 165 0]/255); 
    end
    scatter( intersectvdxy(i,1) , intersectvdxy(i,2), 'filled', 'MarkerFaceColor', 'm');
    hold off;
    xlim([-5 5]); ylim([-5 5]);
    colorbar;

    fig10 = figure(10); clf(fig10);
    contour( x, y, F, [.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
    hold on;
    for n=1:pref.N
        rectangle('Position',[xy(n,1)-r xy(n,2)-r 2*r 2*r],'Curvature', 1, 'EdgeColor', 'none', 'FaceColor', [0.9 0.9 0.9])
    end
    scatter( xy(:,1) , xy(:,2), [], color, 'filled');
    text( xy(:,1), xy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','black');
    plot( xy(convh,1) , xy(convh,2),'m-');
    scatter( points_inhull{i}(:,1) , points_inhull{i}(:,2), [], 'k', 'filled');
    for w=1:length(wells{i})
       plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [255 165 0]/255); 
    end
    hold off;
    xlim([-5 5]); ylim([-5 5]);
    colorbar;

pause(1)
end



F_max = max(F(:));
F_std = F / F_max; % Normalise distribution so values go from 0-1.
distance = reshape(pdist2(xy, [X(:) Y(:)], 'euclidean', 'Smallest', 1), size(X));
distance_max = max(distance(:));
distance_std = distance / distance_max; % Normalise distances so values go from 0-1.

% The scale between distribution density and distance at the different intersection point.
scale = (F_intersectvdxy/F_max)' ./ (xys_hdist / distance_max);
[~, sort_scale] = sort( scale, 'descend' );

for i=sort_scale
    %i=2;
    r =  xys_hdist(i);
    r_std = r / distance_max;
    Fi = F_intersectvdxy(i);
    Fi_std = Fi / F_max;
    F_visible = (F_std/Fi_std > distance_std/r_std);
    F_new = F_std - (F_std - distance_std) * scale(i) ;

    fig12 = figure(12); clf(fig12);
    hold on;
    for n=1:pref.N
        rectangle('Position', [xy(n,1)-r xy(n,2)-r 2*r 2*r], 'Curvature', 1, 'EdgeColor', 'none', 'FaceColor', [0.9 0.9 0.9]);
    end
    imagesc(x, y, F_new, 'AlphaData', F_visible);
    %plot( xy(convh,1) , xy(convh,2),'m-');
    scatter( intersectvdxy(i,1) , intersectvdxy(i,2), 'filled', 'MarkerFaceColor', 'm');
    for w=1:length(wells{i})
       plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [255 165 0]/255); 
    end
    scatter( xy(:,1) , xy(:,2), 50, color, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 2);
    contour(x, y, distance_std, [xys_hdist(i)/max(distance(:)) xys_hdist(i)/max(distance(:))], 'Color', 'black');
    hold off;
    % Set color
    caxis([-1,1]);
    colormap(fig12, 'default')
    map = colormap;
    map2(33:64,:) = map(1:2:end,:);
    map2(32,:) = [0 0 0];
    map2(1:31,:) = repmat([0.7 0.7 0.7], 31, 1);
    colormap(map2);
    % Misc
    set(gca,'ydir', 'normal');
    xlim([-5 5]); ylim([-5 5]);
    colorbar;
    
pause(0.3);
end



% 
% figure(13);
% imagesc(x, y, distance);
% hold on;
% scatter( intersectvdxy(i,1) , intersectvdxy(i,2), 'k*');
% contour(x, y, distance_std, [xys_hdist(i)/max(distance(:)) xys_hdist(i)/max(distance(:))], 'Color', 'black');
% hold off;
% xlim([-5 5]); ylim([-5 5]);
% set(gca,'ydir', 'normal');
% colormap(flipud(winter(25)));
% colorbar;
% 
% 
% figure(14);
% %imagesc(x, y, F);
% imagesc(x, y, F_std);
% colorbar;
% xlim([-5 5]); ylim([-5 5]);
% set(gca,'ydir', 'normal');
% 


