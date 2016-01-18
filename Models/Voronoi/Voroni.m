%% Voronoi diagram
% version 0.01
% Jonas K. Sekamane
%
% *BOUNDED VORONOI DIAGRAM*
% Extending the the _voronoi_ and _voronoin_ functions with a boundary, 
% such that any open regions can be patched.

clearvars;

%% Setup

% Points - Random firm locations
n = 1;
%x0 = rand(1, n)*10-5; % Initial x-position of firm
%y0 = rand(1, n)*10-5; % Initial y-position of firm
%xy = [x0' y0'];
[x0, y0] = pol2cart( rand(n,1)*2*pi , rand(n,1)*3 );
xy = [x0 y0];

if (size(unique(xy,'rows'),1) < n)
    warning('Less unique xy-firm coordinates than number of firms.');
    % TO-DO: Test how often likely this situation is, before providing fix.
end
    
%% Boundary box

% Boundary box
bbox = [-5 5 5 -5 -5; 5 5 -5 -5 5]';

% Calculate the maximum distance between any two points within the boundary
bdistance = sqrt((-5-5)^2+(-5-5)^2)*2; % Twice the length of diagonal for a square boundary box.
% For more general convex boundary use covhull() + algorithm described here: http://stackoverflow.com/a/2736383

% Boundary box line segments
bvx = [ bbox(1:end-1,1)'; circshift( bbox(1:end-1,1)', [0,-1] ) ];
bvy = [ bbox(1:end-1,2)'; circshift( bbox(1:end-1,2)', [0,-1] ) ];


%% Special case (n=1) then area is given by boundary
if(n==1)
    BXY = [bvx' bvy']; % [x1 x2 y1 y2]
    BXY = BXY(:, [1 3 2 4]); % [x1 y1 x2 y2]
    all = reshape(BXY', 2, size(BXY,1)*2)';
    C{1} = 1:length(all);
    
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
        vx_extent = vx;
        vy_extent = vy;
    else
        % Default functions
        [vx,vy] = voronoi(xy(:,1)' , xy(:,2)');
        % Coordinates of the inner Voronoi edges/end points
        %[V1,CB] = voronoin(xy);
    end

    % Default plot
    figure(1);
    % voronoi(xy(:,1)' , xy(:,2)');
    plot(vx,vy, 'b');
    xlim([-6 6]); ylim([-6 6]);
    hold on;
    s = scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
    hold off;



    %% Extend lines
    % Extend the lines to insure that they intersect with the boundary box.
    % Only extend lines where the endspoints are unconnected with other lines 
    % (and extend in the direction away from the connected end). The length 
    % after extension should be the diagonal length of the boundary box.

    if(n>2)

        % Identify the Voronoi edges/endpoints that are within the boundary.
        [~, vui] = unique(vx(:)); 
        V0 = [vx(sort(vui)) vy(sort(vui))];
        in = inpolygon(V0(:,1), V0(:,2), bbox(:,1), bbox(:,2));


        % Identify Voronoi line segments with common endpoints inside boundary.
        vend = NaN(size(vx));
        for line=1:length(vx);
            % If line is connected to any other line (ie. 1:N~=line).
            others = find( 1:length(vx) ~= line );
            for i=1:2
                 connected = ismember(vx(i,line), vx(:,others));
                 outside = ismember(vx(i,line), V0(~in,1));
                 vend(i,line) = connected+outside;
            end
        end
        % Find the Voronoi lines that have exactly one shared endpoint within
        % boundary.
        count = sum(logical(vend),1);
        lines = find(count==1);
        % Loop through all lines with one common endpoint (in boundary) and
        % extend its length (to insure that the line intersects with the boundary)
        vx_extent = vx;
        vy_extent = vy;
        for line = lines
            % The common endpoint
            i = find(vend(:,line));
            % The other endpoint (opposite of i)
            j = 1+abs(2-i); % Takes value 2 if i=1 and takes value 1 if i=2.
            % Extend the line in the direction away from the common endpoint
            [direction, ~] = cart2pol( vx(j,line)-vx(i,line), vy(j,line)-vy(i,line) );
            % Setting the length equal to the maximumum distance between any two
            % points within boundary.
            [dx, dy] = pol2cart(direction, bdistance);
            % New end point
            vx_extent(j,line) = vx(i,line) + dx;
            vy_extent(j,line) = vy(i,line) + dy;
        end

    end

    %% Line intersection

    % Formatting to fit lineSegmentIntersect
    VXY = [vx_extent' vy_extent']; % [x1 x2 y1 y2]
    VXY = VXY(:, [1 3 2 4]); % [x1 y1 x2 y2]
    BXY = [bvx' bvy']; % [x1 x2 y1 y2]
    BXY = BXY(:, [1 3 2 4]); % [x1 y1 x2 y2]

    out = lineSegmentIntersect(VXY, BXY);


    %% Split lines at intersection points

    % New variables that simply later code
    A = out.intAdjacencyMatrix;
    IX = out.intMatrixX;
    IY = out.intMatrixY;
    s = size(A);

    % Create new matrix for the split line segments
    VXY_split = NaN(s(1)+sum(A(:)), 4);
    % Count the number of intersections for each voronoi line.
    vbreaks = sum(A,2);
    % Loop through all voronoi lines that has an intersection.
    v = 1;
    for vline = find(vbreaks)'

        % Find index for intersection points.
        bline = find(A(vline,:));
        % Create matrix with all points. Each row is a point with form [x y].
        % So need to reshape VXY. Each split gives two points; start and end.
        points = [ reshape( VXY(vline,:), 2, 2)'; ...
                   repmat( [IX(vline,bline)' IY(vline,bline)'], 2, 1) ];
        % All points on same line; sort primarily by x and secondarily by y.
        pointssorted = sortrows(points);
        % Reformat so row is line with form [x1 y1 x2 y2]. New length is l.
        l = vbreaks(vline)+1;
        VXY_split(v:v-1+l,:) = reshape(pointssorted', 4, l)';
        v = v+l;

    end
    VXY_split(v:end,:) = VXY(vbreaks==0,:); % Lines with no split

    % Create new matrix for the split line segments
    BXY_split = NaN(s(2)+sum(A(:)), 4);
    % Count the number of intersections for each boundary line.
    bbreaks = sum(A,1);
    % Loop through all boundary lines that has an intersection.
    b = 1;
    for bline = find(bbreaks)

        % Find index for intersection points.
        vline = find(A(:,bline));
        % Create matrix with all points. Each row is a point with form [x y].
        % So need to reshape BXY. Each split gives two points; start and end.
        points = [ reshape( BXY(bline,:), 2, 2)'; ...
                   repmat( [IX(vline,bline) IY(vline,bline)], 2, 1) ];
        % All points on same line; sort primarily by x and secondarily by y.
        if mod(bline,2) % even number
            pointssorted = sortrows(points, [1 2]);
        else
            pointssorted = sortrows(points, [2 1]);
        end
        % Reformat so row is line with form [x1 y1 x2 y2]. New length is l.
        l = bbreaks(bline)+1;
        BXY_split(b:b-1+l,:) = reshape(pointssorted', 4, l)';
        b = b+l;

    end
    BXY_split(b:end,:) = BXY(bbreaks==0,:); % Lines with no split

    % Delete lines outside boundary (.0000000000001 fixes rounding error).
    VXY_split(find( sum( (abs(VXY_split)>5.0000000000001) ,2) ),:) = [];


    %% Match lines with respective firm

    % Calculating distance from line segments to points

    % Distance to mid-point of line segment
    BXY_split_mid = [mean(BXY_split(:,[1 3]),2) mean(BXY_split(:,[2 4]),2)] ;
    bdist_mid = pdist2(xy, BXY_split_mid);
    [~, IB] = min(bdist_mid);

    % bdist = pldist2(xy, BXY_split);
    % % Find closest firm for each line segment
    % [~, IB] = min(bdist);

    VXY_split_mid = [mean(VXY_split(:,[1 3]),2) mean(VXY_split(:,[2 4]),2)] ;
    vdist_mid = pdist2(xy, VXY_split_mid);
    [~, IV] = sort(vdist_mid, 1);
    IV(3:end,:) = [];

    % vdist = pldist2(xy, VXY_split);
    % % Each voronoi line is borders two firms. Keep the two with min distance.
    % [~, IV] = sort(vdist, 1);
    % IV(3:end,:) = [];

    % Map the line segment to the firms.
    all_split = [BXY_split; VXY_split];
    C_all = num2cell(repmat(NaN,1,n)');
    for i=1:length(C_all)
        bline = find(IB==i);
       [~, vline] = find((IV==i));
       C_all{i} = [ bline vline'+size(BXY_split,1) ];
    end

    %% Reformat from lines to points 
    all = reshape(all_split', 2, size(all_split,1)*2)';
    C = cellfun(@(x) [x*2-1 x*2], C_all, 'UniformOutput', false);

end

%% Convex hull / Polygon

% Arrange points counterclock wise around convex hull
hull = num2cell(repmat(NaN,1,n)');
for i = 1:n
    hull_i = convhull( all(C{i},1), all(C{i},2) );
    hull{i} = hull_i;
end


%% Figures
figure(2);
plot( VXY_split(:,[1 3])', VXY_split(:,[2 4])');
xlim([-6 6]); ylim([-6 6]);
hold on;
plot( BXY_split(:,[1 3])', BXY_split(:,[2 4])');
text( BXY_split(:,1)'+0.1, BXY_split(:,2)'+0.1, cellstr(num2str([1:length(BXY_split)]')) );
scatter(IX(A),IY(A),[],'k','filled');
%scatter(V1(:,1), V1(:,2));
%scatter(bbox(1:end-1,1),bbox(1:end-1,2));
%text(bbox(1:end-1,1)-0.1, bbox(1:end-1,2)-0.1, cellstr(num2str([1:4]')) );
scatter( xy(:,1)' , xy(:,2)', 'filled');
text(xy(:,1)'+0.1, xy(:,2)'+0.1, cellstr(num2str([1:n]')) );
%plot(linex,liney,':');
%scatter(xi_temp, yi_temp);
%patch(test(:,1),test(:,2),3);
hold off;
% xlabel('x'); ylabel('y');
% %voronoi(xy(:,1)' , xy(:,2)');
% hold off;


% figure(3);
% plot(vx,vy);
% xlim([-10 10]); ylim([-10 10]);
% hold on;
% plot(bbox(:,1),bbox(:,2));
% hold off;
% 
% figure(31);
% plot(vx_extent,vy_extent);
% xlim([-10 10]); ylim([-10 10]);
% hold on;
% plot(bbox(:,1),bbox(:,2));
% hold off;

% figure(4);
% plot(vx(:,bbii(1,2)),vy(:,bbii(1,2)));
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% plot( bbox(bbii(1,1),1) , bbox(bbii(1,1),2) );
% hold off;
% bbii(1,1)
% 

% figure(33);
% plot(vx_extent,vy_extent);
% xlim([-8 8]); ylim([-8 8]);
% hold on;
% plot(bbox(:,1),bbox(:,2));
% scatter(IX(A),IY(A),[],'r');
% hold off;

% figure(5);
% plot( BXY_split(:,[1 3])', BXY_split(:,[2 4])');
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% plot( VXY_split(:,[1 3])', VXY_split(:,[2 4])');
% hold off;

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


figure(555);
scatter( xy(:,1)' , xy(:,2)', 'filled');
hold on;
if(n>1) plot( all_split(:,[1 3])', all_split(:,[2 4])', 'k'); end
%patch( all_split(m,[1 3])', all_split(m,[2 4])', 'red');
for i = 1:n
    X = all(C{i},1);
    Y = all(C{i},2);
    patch(X(hull{i}), Y(hull{i}), i, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
xlim([-6 6]); ylim([-6 6]);
hold off;
