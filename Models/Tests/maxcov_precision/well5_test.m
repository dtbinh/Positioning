function [xy_centroid_well, xy_well_idx] = well5_test( firms, customers, F )

%X = customers(:,1);
%Y = customers(:,1);
xy = firms;
pref.N = size(xy,1);

% TO-DO: Create population function with input [x y pref] and output F or [F, F_l, F_r]

% Boundary box
bbox = [-5 5 5 -5 -5; 5 5 -5 -5 5]';
bboxxy = [bbox(1:end-1,:) circshift(bbox(1:end-1,:),-1)];

% %% Population distribution
% pref.boundary = 10; % Number of standard deviations
% pref.resolution = 50; % Length of the square (even number to include (0,0))
% b = pref.boundary/2;
% [x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
% [X,Y] = meshgrid(x,y);
% %Population
% pref.mu = 1.5; % Mean of subpopulation
% pref.n_ratio = 2; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
% sd = 0.5; % Standard deviation of each subpopulation
% mu_r = [pref.mu 0]; % Subpopulation mean only deviate on x-axis
% sigma_r = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
% F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r); % Subpopulation pdf evaluated at each point in grid/square
% F_r = reshape(F_r,length(y),length(x)); % Formated as grid
% mu_l = [-pref.mu 0]; % Subpopulation mean only deviate on x-axis
% sigma_l = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
% F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l); % Subpopulation pdf evaluated at each point in grid/square
% F_l = reshape(F_l,length(y),length(x)); % Formated as grid
% weight = pref.n_ratio/(1+pref.n_ratio);
% F = (F_l + F_r*1/pref.n_ratio)/4;


%% Loation of firms
% 
% pref.N = 8;
% [x0, y0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*3 );
% xy = [x0 y0];
% color = brewermap(pref.N,'Paired');


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
%F_v3xy = (mvnpdf(vxy, mu_l, sigma_l) + mvnpdf(vxy, mu_r, sigma_r)*1/pref.n_ratio)/4;


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
%F_intersectvdxy = (mvnpdf(intersectvdxy, mu_l, sigma_l) + mvnpdf(intersectvdxy, mu_r, sigma_r)*1/pref.n_ratio)/4;


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


% Boundary solutions (outside convex hull)


% A vector for each face of the polygon
convh_vector = diff(xy(convh,:)); % Arranged CCW
% Convex hull points including mid-points (half vector)
%convhmid = [xy(convh,:); xy(convh(1:end-1), :) + convh_vector/2];
convhmid = [xy(convh(1:end-1), :) + convh_vector/2];
% Find the point on the convex hull that is furthest from the boundary (closest to zero)
[convhdist, convhmin] = pdist2(convhmid, [0 0], 'euclidean', 'Smallest', 1);

% Is zero part of the convex hull.
%contained = InPolygon(0, 0, xy(convh,1) , xy(convh,2));
%sign = contained*2-1;
%sign = 1;
%nch = size(xy(convh,:),1);
nch=0;
%smalloffset = 0.01; % 0.1 std. dev.

% if (convhmin <= nch)
%     % Test if it is a corner point on the convex hull 
%     
%     % New point at
%     new_xy = convhmid(convhmin,:) * (1+smalloffset*sign);% ;
%     
%     % Vector orthogonal to the line to zero, ie. new face of voronoi region
%     parallel_vector = convhmid(convhmin,:) * [0 -1; 1 0];%([0 -1; 1 0]*contained + [1 0; 0 -1]*(1-contained));
%     parallel_norm = sqrt(dot(parallel_vector, parallel_vector, 2));
%     polypoints = [new_xy - parallel_vector/parallel_norm*sqrt(200)*2; new_xy + parallel_vector/parallel_norm*sqrt(200)*2];
%     ployxy = reshape(polypoints', 4, [])';
%     
%     % Extend the voronoi lines to insure that they intersect the boundary.
%     [vx_extend, vy_extend] = voronoi_extend(vx, vy, bbox, sqrt(200)*2);
%     
%     % Find the voronoi lines/region surrounding point
%     VXY_mid = [mean(vx_extend)' mean(vy_extend)'];
%     [vdist_mid, vpoint_mid]  = pdist2(xy, VXY_mid, 'euclidean', 'Smallest', 2);
%     [~, vline] = find(vpoint_mid == convh(convhmin));
%     pointvline = [vx_extend(1,vline)' vy_extend(1,vline)' vx_extend(2,vline)' vy_extend(2,vline)'];
%     
%     bboxintersect = lineSegmentIntersect(pointvline, bboxxy);
%     bboxintersectxy = [bboxintersect.intMatrixX(bboxintersect.intAdjacencyMatrix) bboxintersect.intMatrixY(bboxintersect.intAdjacencyMatrix)];
%     
%     % Find intersect points between the existing Voronoi lines and the new face of voronoi region and the boundary points.
%     %polyintersect = lineSegmentIntersect([ployxy; bboxxy], pointvline);
%     %polyintersectxy = [polyintersect.intMatrixX(polyintersect.intAdjacencyMatrix) polyintersect.intMatrixY(polyintersect.intAdjacencyMatrix)];
%     polyintersect = lineSegmentIntersect(ployxy, bboxxy);
%     polyintersectxy = [polyintersect.intMatrixX(polyintersect.intAdjacencyMatrix)' polyintersect.intMatrixY(polyintersect.intAdjacencyMatrix)'];
%     
%     polyvintersect = lineSegmentIntersect(ployxy, pointvline);
%     polyvlineintersectxy = [polyvintersect.intMatrixX(polyvintersect.intAdjacencyMatrix)' polyvintersect.intMatrixY(polyvintersect.intAdjacencyMatrix)'];
%     
%     %vintersect = any(polyintersect.intAdjacencyMatrix(2:end,:));
%     %pointvline(~vintersect,:)
%     %pointvline(vintersect,:)
%     vintersectlines = reshape( pointvline', 2, [])';
%     [vlineintersectlinesu,~, lineic] = unique(vintersectlines, 'rows');
%     nlines = histc(lineic, unique(lineic));
%     vlineintersectlinesxy = vlineintersectlinesu(find(nlines >= 2),:);
%     
% %     potentialxy = [bbox(1:4,:); bboxintersectxy; polyintersectxy; vlineintersectlinesxy];
% %     %[new_xy; xy]
% %     [newcornersdist, newcorners_all] = pdist2([new_xy; xy], potentialxy, 'euclidean', 'Smallest', 2);
% %     % Find all the potential points where the new points is closest (or second closest)
% %     newcorners = any(newcorners_all == 1);
% %     % Adjusting the bbox points.
% %     unequaldist = find(diff(newcornersdist(:,1:4), [], 1) > 0); % If equal distance then keep potential point, otherwise ...
% %     % Otherwise only keep the bbox point if it is the closest.
% %     newcorners(unequaldist) = newcorners(unequaldist).*(newcorners_all(1,unequaldist) == 1);
% 
%     potentialxy = [bboxintersectxy; polyintersectxy; vlineintersectlinesxy];
%     %[new_xy; xy]
%     [newcornersdist, newcorners_all] = pdist2([new_xy; xy], potentialxy, 'euclidean', 'Smallest', 2);
%     % Find all the potential points where the new points is closest (or second closest)
%     newcorners = any(newcorners_all == 1);
%     % Adjusting the bbox points.
%     unequaldist = find(diff(newcornersdist, [], 1) > 0); % If equal distance then keep potential point, otherwise ...
%     % Otherwise only keep the bbox point if it is the closest.
%     newcorners(unequaldist) = newcorners(unequaldist).*(newcorners_all(1,unequaldist) == 1);
%     potentialxy_temp = potentialxy(newcorners,:);
%     newbbox = any(abs(potentialxy_temp) == 5,2);
%     
%     [~, bboxcorner] = pdist2(bbox(1:4,:), potentialxy_temp(newbbox,:), 'euclidean', 'Smallest', 2);
%     bboxcorneru = unique(bboxcorner);
%     bboxn = histc(bboxcorner(:), unique(bboxcorner));
%     if diff(bboxcorner(1,:)) == 0
%        % If they share closest corner
%        bboxcorner2 = unique(bboxcorner(1,:));
%     elseif length(bboxcorneru) == 2
%         % If they share both closest corners, then no bbox corners needed.
%         bboxcorner2 = [];
%     elseif any(bboxn-1)
%        % If they share corner
%        bboxcorner2 = bboxcorneru(bboxn == 2,:);
%     else
%        % Otherwise pick closest corners.
%        bboxcorner2 = bboxcorner(1,:);
%     end
%     
%     %convh(convhmin)
%     
%     %pointvline(any(polyintersect.intAdjacencyMatrix),:)
%     % inserct intersection points, and include all vlines around point
%     % pointvline(any(polyintersect.intAdjacencyMatrix),:)
%     
%     % Keep the points where any of the points voronoi lines intersect or where the points of boundary intersect.
%     %vvbb = [pointvline(any(polyintersect.intAdjacencyMatrix),:); bboxxy(any(polyintersect.intAdjacencyMatrix(2:end,:), 2),:)];
%     %vlineintersect = reshape( vvbb', 2, [])';
%     %[vlineintersectu,~, ic] = unique(vlineintersect, 'rows');
%     %n = histc(ic, unique(ic));
%     %vlineintersectxy = vlineintersectu(find(n >= 2),:);
%     %if ~isempty(vlineintersectlinesxy)
%     %   [~, indx] = ismember(vlineintersectxy, vlineintersectlinesxy, 'rows');
%     %   vlineintersectxy(find(indx), :) = [];
%     %end
%     potentialxy_all = potentialxy(newcorners,:);
%     inside = ~any(abs(potentialxy_all)>5+0.001,2);
%     
%     % Combine points and order CCW as polygon/convex hull
%     polyraw = [potentialxy_all(inside,:); polyvlineintersectxy; bbox(bboxcorner2,:)];
%     polyorder = convhull(polyraw(:,1), polyraw(:,2));
%     poly = polyraw(polyorder,:);
%     
%     
% else
    % Else it is a mid-point on the convex hull
    
    % New point at mid-point plus vector rotated 90 degrees CCW.
    new_vector = convh_vector(convhmin-nch, :) * [0 -1; 1 0];
    new_vector_norm = sqrt(dot(new_vector, new_vector, 2)); % Norm of vector
    new_xy = convhmid(convhmin,:) + new_vector/new_vector_norm * convhdist;
    
    % The two new faces of the polygon at +/- 45 degree angle of the convex hull
    ccw = pi/4;
    cw = 2*pi-ccw;
    poly_ccw = convh_vector(convhmin-nch, :) * [cos(ccw) -sin(ccw); sin(ccw) cos(ccw)];
    poly_ccw_norm = sqrt(dot(poly_ccw, poly_ccw, 2));
    poly_cw = convh_vector(convhmin-nch, :) * [-cos(cw) sin(cw); -sin(cw) -cos(cw)];
    %poly_cw = poly_ccw * [0 -1; 1 0];
    poly_cw_norm = sqrt(dot(poly_cw, poly_cw, 2));
    polypoints = repmat(convhmid(convhmin,:),2) + [poly_ccw/poly_ccw_norm 0 0; 0 0 poly_cw/poly_cw_norm]*sqrt(200)*2; %
    
    ployxy = reshape(polypoints', 4, [])';
    % Find the intersection points with the boundary
    polyintersect = lineSegmentIntersect(ployxy, bboxxy);
    polyintersectxy = [polyintersect.intMatrixX(polyintersect.intAdjacencyMatrix) polyintersect.intMatrixY(polyintersect.intAdjacencyMatrix)];
    
    % Decide which bbox corners to include (1, 2, or none).
    [~, bboxcorner] = pdist2(bbox(1:4,:), polyintersectxy, 'euclidean', 'Smallest', 2);
    bboxcorneru = unique(bboxcorner);
    bboxn = histc(bboxcorner(:), unique(bboxcorner));
    if diff(bboxcorner(1,:)) == 0
       % If they share closest corner
       bboxcorner2 = unique(bboxcorner(1,:));
    elseif length(bboxcorneru) == 2
        % If they share both closest corners, then no bbox corners needed.
        bboxcorner2 = [];
    elseif any(bboxn-1)
       % If they share corner
       bboxcorner2 = bboxcorneru(bboxn == 2,:);
    else
       % Otherwise pick closest corners.
       bboxcorner2 = bboxcorner(1,:);
    end
    
    % Combine points and order CCW as polygon/convex hull
    polyraw = [convhmid(convhmin,:); bbox(bboxcorner2,:); polyintersectxy];
    polyorder = convhull(polyraw(:,1), polyraw(:,2));
    poly = polyraw(polyorder,:);
% end

wells{length(wells)+1}{1} = poly;

% Delete empty wells.
wells = wells(~cellfun('isempty',wells));

% Calculate market share for each well
well_centroid = nan(1,2);
well_id = nan(1,2);
well_count = 1;
for i=1:length(wells)
    for w=1:length(wells{i})
       idx = inpolygon(customers(:,1), customers(:,2), wells{i}{w}(:,1), wells{i}{w}(:,2) );
       F_well(well_count) = sum( F(idx) )  ;
       well_centroid(well_count,:) =  ([customers(idx,1) customers(idx,2)]' * F(idx))' ./ F_well(well_count);
       well_id(well_count,:) = [i w];
       well_count = well_count+1;
    end
end

% Remove wells with zero market
wellpositive = (F_well > 0);
F_well(~wellpositive) = [];
well_centroid(~wellpositive,:) = [];

[~, well_max] = max(F_well);
%xy_well_max = well_centroid(well_max,:);

% Function output
xy_centroid_well = well_centroid;
xy_well_idx = well_max;


% %% PLOTS
% 
% 
% figure(11);
% scatter( xy(:,1) , xy(:,2), [], 'b', 'filled');
% hold on;
% plot(vx,vy, 'b');
% %plot(vx_extend,vy_extend, 'b');
% %triplot( DT, ':', 'Color', [0.5 0.5 0.5]);
% scatter( 0, 0, [], 'r*');
% nr = length(wells);
% j = nr;
% for i=1:nr    
%     for w=1:length(wells{i})
%        plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [1 0.65 (1-j/(nr+1))]); 
%     end
%     j = j-1;
% end
% %plot( poly(:,1) , poly(:,2), 'Color', [1 0.65 0]);
% %text( xy(:,1), xy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','black');
% scatter( new_xy(1) , new_xy(2), [], 'k', 'filled');
% %scatter( polyraw(:,1), polyraw(:,2));
% %scatter( potentialxy(:,1), potentialxy(:,2), 'r');
% %plot( polypoints(:,1) , polypoints(:,2), 'k');%, 'Color', [1 0.65 0]);
% hold off;
% xlim([-5 5]); ylim([-5 5]);
% colorbar;
% 
% 
% figure(15);
% scatter( 0, 0, [], 'r*');
% hold on;
% nr = length(sort_dist);
% j = length(sort_dist);
% for i=sort_dist    
%     for w=1:length(wells{i})
%        plot( wells{i}{w}(:,1) , wells{i}{w}(:,2), 'Color', [1 0.65 (1-j/(nr+1))]); 
%     end
%     j = j-1;
% end
% scatter( convhmid(:,1) , convhmid(:,2), [], 'm', 'filled');
% scatter( convhmid(convhmin,1) , convhmid(convhmin,2), [], 'k');
% scatter( new_xy(1) , new_xy(2), [], 'k', 'filled');
% plot( xy(convh,1) , xy(convh,2),'m-');
% plot( poly(:,1) , poly(:,2), 'Color', [1 0.65 0]);
% %scatter( polyraw(:,1), polyraw(:,2));
% hold off;
% xlim([-5 5]); ylim([-5 5]);
% colorbar;
% 
% 
end