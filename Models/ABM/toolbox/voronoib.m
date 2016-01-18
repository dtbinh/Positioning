function [P, C, H] = voronoib(xy)
%VORONOIB
%   An extension of matlab's _voronoi_ function to include a boundary. Ie.
%   open regions are patched and works with both 1 and 2 generating points.
%   
%   Jonas K. Sekamane. 
%   Version 0.01
%

%clearvars;

%% Setup

% Points - Random firm locations
%n = 3;
%x0 = rand(1, n)*10-5; % Initial x-position of firm
%y0 = rand(1, n)*10-5; % Initial y-position of firm
%xy = [x0' y0'];
%[x0, y0] = pol2cart( rand(n,1)*2*pi , rand(n,1)*3 );
%xy = [x0 y0];
n = size(xy,1);

if (size(unique(xy,'rows'),1) < n)
    warning('Less unique xy-firm coordinates than number of firms.');
    % TO-DO: Test how often likely this situation is, before providing fix.
end
    
%% Boundary box

% Boundary box
%bbox = [x([1 end end 1 1]); y([end end 1 1 end])]';
bbox = [-5 5 5 -5 -5; 5 5 -5 -5 5]';

% Boundary box line segments
bvx = [ bbox(1:end-1,1)'; circshift( bbox(1:end-1,1)', [0,-1] ) ];
bvy = [ bbox(1:end-1,2)'; circshift( bbox(1:end-1,2)', [0,-1] ) ];

% Calculate the maximum distance between any two points within the boundary
%bdistance = sqrt((pref.boundary)^2+(pref.boundary)^2)*2; % Twice the length of diagonal for a square boundary box.
bdistance = sqrt((-5-5)^2+(-5-5)^2)*2; % Twice the length of diagonal for a square boundary box.
% For more general convex boundary use covhull() + algorithm described here: http://stackoverflow.com/a/2736383


%% Special case (n=1) then market polygon is given by boundary
if(n==1)
    BXY = [bvx' bvy']; % [x1 x2 y1 y2]
    BXY = BXY(:, [1 3 2 4]); % [x1 y1 x2 y2]
    all = reshape(BXY', 2, size(BXY,1)*2)';
    C_temp{1} = 1:length(all);
    
else 
    
    %% Create Voronoi lines
    if(n==2)
        % Normal vector between the two points (going though midpoint). First
        % term is normal with length 1, second term increases the distance to 
        % bdistance, insuring intersection with boundary.
        xyn = null(diff(xy))' .* bdistance/2;
        xynormal = [mean(xy)-xyn; mean(xy)+xyn];
        vx = xynormal(:,1);
        vy = xynormal(:,2);
        vx_extend = vx;
        vy_extend = vy;
    else
        % Default functions
        [vx,vy] = voronoi(xy(:,1)' , xy(:,2)');
        % Coordinates of the inner Voronoi edges/end points
        %[V1,CB] = voronoin(xy);
    end

%     % Default plot
%     figure(1);
%     % voronoi(xy(:,1)' , xy(:,2)');
%     plot(vx,vy, 'b');
%     xlim([-6 6]); ylim([-6 6]);
%     hold on;
%     s = scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
%     hold off;



    %% Extend lines
    % Extend the lines to insure that they intersect with the boundary box.
    % Only extend lines where the endspoints are unconnected with other lines 
    % (and extend in the direction away from the connected end). The length 
    % after extension should be the diagonal length of the boundary box.

    if(n>2)

        [vx_extend, vy_extend] = voronoi_extend(vx, vy, bbox, bdistance);

    end

    %% Line intersection

    % Formatting to fit lineSegmentIntersect
    VXY = [vx_extend' vy_extend']; % [x1 x2 y1 y2]
    VXY = VXY(:, [1 3 2 4]); % [x1 y1 x2 y2]
    BXY = [bvx' bvy']; % [x1 x2 y1 y2]
    BXY = BXY(:, [1 3 2 4]); % [x1 y1 x2 y2]

    out = lineSegmentIntersect(VXY, BXY);


    %% Split lines at intersection points

    [VXY_split, BXY_split] = voronoi_split(VXY, BXY, out.intAdjacencyMatrix, out.intMatrixX, out.intMatrixY);

    % Delete lines outside boundary (.0000000000001 fixes rounding error).
    VXY_split(find( sum( (abs(VXY_split) > abs(BXY(1,1))+.0000000000001) ,2) ),:) = [];


    %% Match lines with respective firm

    % Calculating distance from line segments to points

    % Distance to mid-point of line segment
    BXY_split_mid = [ mean( BXY_split(:,[1 3]), 2) mean( BXY_split(:,[2 4]), 2) ] ;
    bdist_mid = pdist2(xy, BXY_split_mid);
    [~, IB] = min(bdist_mid);

    VXY_split_mid = [ mean( VXY_split(:,[1 3]), 2) mean( VXY_split(:,[2 4]), 2) ] ;
    vdist_mid = pdist2(xy, VXY_split_mid);
    % Each voronoi line is the border between two firms. Keep the two with min distance.
    [~, IV] = sort(vdist_mid, 1);
    IV(3:end,:) = [];

    % Map the line segment to the firms.
    all_split = [BXY_split; VXY_split];
    C_all = num2cell(NaN(n,1));
    for i=1:length(C_all)
        bline = find(IB==i);
       [~, vline] = find((IV==i));
       C_all{i} = [ bline vline'+size(BXY_split,1) ];
    end

    %% Reformat from lines to points 
    all = reshape(all_split', 2, size(all_split,1)*2)';
    C_temp = cellfun(@(x) [x*2-1 x*2], C_all, 'UniformOutput', false);

end

%% Convex hull / Polygon

% Arrange points counterclock wise around convex hull
hull = num2cell(NaN(n,1));
for i = 1:n
    hull_i = convhull( all(C_temp{i},1), all(C_temp{i},2) );
    hull{i} = hull_i;
end



%% Output variable
P = all; % Points
C = C_temp; % Connecting each row in xy to the points in all.
H = hull; % Convex hulls


end


%% Figures
% figure(2);
% plot( VXY_split(:,[1 3])', VXY_split(:,[2 4])');
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% plot( BXY_split(:,[1 3])', BXY_split(:,[2 4])');
% text( BXY_split(:,1)'+0.1, BXY_split(:,2)'+0.1, cellstr(num2str([1:length(BXY_split)]')) );
% scatter(out.intMatrixX(out.intAdjacencyMatrix),out.intMatrixY(out.intAdjacencyMatrix),[],'k','filled');
% scatter( xy(:,1)' , xy(:,2)', 'filled');
% text(xy(:,1)'+0.1, xy(:,2)'+0.1, cellstr(num2str([1:n]')) );
% hold off;
%
% figure(3);
% plot(vx,vy);
% xlim([-10 10]); ylim([-10 10]);
% hold on;
% plot(bbox(:,1),bbox(:,2));
% hold off;
% 
% figure(31);
% plot(vx_extend,vy_extend);
% xlim([-10 10]); ylim([-10 10]);
% hold on;
% plot(bbox(:,1),bbox(:,2));
% hold off;
%
% figure(4);
% plot(vx(:,bbii(1,2)),vy(:,bbii(1,2)));
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% plot( bbox(bbii(1,1),1) , bbox(bbii(1,1),2) );
% hold off;
% bbii(1,1)
% 
% figure(33);
% plot(vx_extend,vy_extend);
% xlim([-8 8]); ylim([-8 8]);
% hold on;
% plot(bbox(:,1),bbox(:,2));
% scatter(out.intMatrixX(out.intAdjacencyMatrix),out.intMatrixY(out.intAdjacencyMatrix),[],'r');
% hold off;
%
% figure(5);
% plot( BXY_split(:,[1 3])', BXY_split(:,[2 4])');
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% plot( VXY_split(:,[1 3])', VXY_split(:,[2 4])');
% hold off;
%
% %m = C{4}; 
% figure(55);
% scatter( xy(:,1)' , xy(:,2)', 'filled');
% hold on;
% plot( all_split(:,[1 3])', all_split(:,[2 4])', 'k');
% %patch( all_split(m,[1 3])', all_split(m,[2 4])', 'red');
% for i = 1:length(C_all)
%     X = all_split(C_all{i},[1 3])';
%     X2 = X(:);
%     Y = all_split(C_all{i},[2 4])';
%     Y2 = Y(:);
%     [X3, Y3] = poly2ccw(X2, Y2);
%     hull = convhull(X2,Y2);
%     %[lat, lon] = polymerge(X(:), Y(:))
%     %patch(all_split(C{i},[1 3])', all_split(C{i},[2 4])', i); % use color i.
%     %patch(lat, lon, i); % use color i.
%     patch(X2(hull), Y2(hull), i, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % use color i.
%     %patch(X3, Y3, i, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % use color i.
% end
% xlim([-6 6]); ylim([-6 6]);
% hold off;
%
% figure(555);
% scatter( xy(:,1)' , xy(:,2)', 'filled');
% hold on;
% if(n>1) plot( all_split(:,[1 3])', all_split(:,[2 4])', 'k'); end
% %patch( all_split(m,[1 3])', all_split(m,[2 4])', 'red');
% for i = 1:n
%     X = all(C{i},1);
%     Y = all(C{i},2);
%     patch(X(hull{i}), Y(hull{i}), i, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% end
% xlim([-6 6]); ylim([-6 6]);
% hold off;
