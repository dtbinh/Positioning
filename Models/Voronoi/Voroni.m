%% Voronoi diagram
% version 0.01
% Jonas K. Sekamane
%
% *BOUNDED VORONOI DIAGRAM*
% Extending the the _voronoi_ and _voronoin_ functions with a boundary, 
% such that any open regions can be patched.

clearvars;

% Points - Random firm locations
n = 3;
%x0 = rand(1, n)*10-5; % Initial x-position of firm
%y0 = rand(1, n)*10-5; % Initial y-position of firm
%xy = [x0' y0'];
[x0, y0] = pol2cart( rand(n,1)*2*pi , rand(n,1)*3 );
xy = [x0 y0];

% Default plot
figure(1);
voronoi(xy(:,1)' , xy(:,2)');
xlim([-6 6]); ylim([-6 6]);

% Default functions
[vx,vy] = voronoi(xy(:,1)' , xy(:,2)');
% Coordinates of the inner Voronoi edges/end points
[V1,C] = voronoin(xy);

% Boundary box
bbox = [-5 5 5 -5 -5; 5 5 -5 -5 5]';

% Calculate the maximum distance between any two points within the boundary
bdistance = sqrt((-5-5)^2+(-5-5)^2); % Length of diagonal for a square boundary box.
% For more general convex boundary use covhull() + algorithm described here: http://stackoverflow.com/a/2736383

% Boundary box line segments
bvx = [ bbox(1:end-1,1)'; circshift( bbox(1:end-1,1)', [0,-1] ) ];
bvy = [ bbox(1:end-1,2)'; circshift( bbox(1:end-1,2)', [0,-1] ) ];

% Coordinates of the boundary box cornors
V2 = [bbox(1:end-1,1) bbox(1:end-1,2)];

% For each point calculate the euclidean distance to boundary points.
distance = pdist2(xy, [bbox(1:end-1,1) bbox(1:end-1,2)], 'euclidean');
% Find for each boundary point the minimum distiance (minVal) and the boundary point (minInd)
[minVal, minInd] = min(distance);
% Mapping the boundary point to each firm
% Create emptry cell array where row is point, and array is the corresponding bounardy points in V2
C2 = num2cell(repmat(NaN,1,n)');
for i=1:length(C2)
   C2{i} = find(minInd==i);
end

%% TO-DO:
% Extend the lines to insure that they intersect with the boundary box.
% Only extend lines where both (or just one) ends are unconnected with other lines 
% (and extend in the direction away from the connected end).
% The minimum length of the extended line should be the diagonal length of the boundary box.
% If boundary polygon instead of boundary box, then something similar to regionprops's MajorAxisLength property.

% Identify the Voronoi edges/endpoints that are within the boundary.
%[vxu vxui] = unique(vx(:)); 
%[vyu vyui] = unique(vy(:)); 
%[sort(vxui) sort(vyui)]
%vx(sort(vxui))
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


% Remove lines that are completely outside boundary
%vx_extent(:,find(~count)) = [];
%vy_extent(:,find(~count)) = [];

%% Line intersection

%%% First appoach:

% Line intersect with boundary box
bbx = []; bby = []; bbii = [];
for j = 1:length(vx_extent)
    [xi, yi, ii] = polyxpoly(bbox(:,1), bbox(:,2), vx_extent(:,j), vy_extent(:,j)); % [xi, yi, ii] = 
    bbx = [bbx; xi];
    bby = [bby; yi];
    bbii = [bbii; ii];
end
V3 = [bbx bby]; % Coordinates of the voronoi diagram lines intersection with the boundary box.

%%% Second appoach:

% Find the intersection points between the Voronoi lines and the boundary.
[xi, yi, ii] = polyxpoly(vx_extent, vy_extent, bvx, bvy);
% Create new matrix with boundary line segments that are split at the 
% intersection points. Size of new matrix.
sb = [2 size(ii,1)+size(bvx,2)];
bvx_split = NaN(sb);
bvy_split = NaN(sb);
% Create new matrix with Voroni line segments that are split at the 
% intersection points. Size of new matrix.
sv = [2 size(ii,1)+size(vx_extent,2)];
vx_extent_split = NaN(sv);
vy_extent_split = NaN(sv);
[~, iiv] = sort(ii(:,1));
for i=1:size(ii,1)
    % Select the boundary split points in new matrix.
    % Using index from polyxpoly and scales to fit new matrix dimension
    % (each split adds two endpoints)
    bline = ii(i,2) + (i-1)*2;
    bvx_split(bline+1) = xi(i); % Split x end point
    bvx_split(bline+2) = xi(i); % Split x start point
    bvy_split(bline+1) = yi(i); % Split y end point
    bvy_split(bline+2) = yi(i); % Split y start point
    % Select the Voroni split points in new matrix.
    % Using index from polyxpoly and scales to fit new matrix dimension
    % (each split adds two endpoints)
    vline = ii(iiv(i),1) + (i-1)*2;
    vx_extent_split(vline) = xi(iiv(i)); % Split x end point
    vx_extent_split(vline+1) = xi(iiv(i)); % Split x start point
end
%bvx_split
%bvy_split
% Now that all splits has been added the remaning empty cells are the
% original boundary points.
bvx_split_nan = find(isnan(bvx_split));
bvx_split(bvx_split_nan) = bvx;
bvy_split(bvx_split_nan) = bvy;
vx_extent_split_nan = find(isnan(vx_extent_split));
vx_extent_split(vx_extent_split_nan) = vx_extent;
%bvx_split
%bvy_split

%bvx_new = bvx(ii(:,2)) 
%bvx(ii(:,2)+1)

%[ii(:,2) ii(:,2)+1]
[ii(:,1) ii(:,1)+1]
vx_extent_split
% [sort(ii(:,1)) sort(ii(:,1))+1]
% 
% bvx_split = bvx
% bvx_split(7) = NaN
% 
% [xi, yi, jj] = polyxpoly(bvx(:), bvy(:), vx_extent(:), vy_extent(:));



%%% Third appoach:

% Formatting to fit lineSegmentIntersect
VXY = [vx_extent(1,:)' vy_extent(1,:)' vx_extent(2,:)' vy_extent(2,:)'];
BXY = [bvx(1,:)' bvy(1,:)' bvx(2,:)' bvy(2,:)'];

out = lineSegmentIntersect(VXY, BXY);

sb = [size(BXY,2)+sum(out.intAdjacencyMatrix(:)), 2];
BX_split = NaN(sb);
BX_split(out.intAdjacencyMatrix) = out.intMatrixX(out.intAdjacencyMatrix)
BX_split(circshift(out.intAdjacencyMatrix,[0,1])) = out.intMatrixX(out.intAdjacencyMatrix)
%bshift = [zeros(1,size(out.intAdjacencyMatrix,2)); out.intAdjacencyMatrix];
%BX_split(find(bshift)) = out.intMatrixX(out.intAdjacencyMatrix)
BX_split(~out.intAdjacencyMatrix) = BXY(~out.intAdjacencyMatrix)

% Loop through all boundary lines
% bintersets = sum(out.intAdjacencyMatrix,2);
% for bi=1:size(BXY,1)
%     if(~bintersets(bi))
%         % No intersects points
%         BXY_split = [BXY_split; BXY(bi,:)]
%     else
%        % intersect points
%        
%     end
%     
% end

%out.intAdjacencyMatrix

%% 

% For each point calculate the euclidean distance to boundary points.
distance2 = pdist2(xy, [bbx bby], 'euclidean');
[B,I] = sort(distance2,1); % find the smallest distances
I = I(1:2,:); % each point has two boundary points. 
%% TO-DO: 
% This should not find 2 points, but find all cases where the distance is equal.
% Loop though and find all points where distance equals minimum distance B.
% SECOND APPROACH; Only look at points in an open region: if any(C{i}~=1) 
% Create emptry cell array where row is point, and array is the corresponding bounardy points in V3
C3 = num2cell(repmat(NaN,1,n)');
for i=1:length(C3)
   C3{i} = [find(I(1,:)==i) find(I(2,:)==i)]; % the two boundary points
end

%% Combine the points:
% TO-DO:
polygon = [V1; V2; V3];

% First appoach:
% polygon2 = polygon(3:end,:);
% C_final2{1} = [1 3 5 10 11];
% C_final2{2} = [1 2 3 8 9];
% C_final2{3} = [2 4 7 8 12]; 
% C_final2{4} = [1 6 9 11];
% C_final2{5} = [2 3 10 12];

% for i = 3:length(C_final2)
%     x = polygon2(C_final2{i},1);
%     y = polygon2(C_final2{i},2);
%     cx = mean(x);
%     cy = mean(y);
%     a = atan2(x - cy, y - cx);
%     [~, order] = sort(a);
%     x = x(order);
%     y = y(order);
%     patch(x,y,i);
% end

% Second appoach:
% for i = 1:length(C)
%     C_final{i} = [C{i} C2{i}+length(V1) C3{i}+length(V1)+length(V2)];
%     
% %     if all(C{i}~=1)   % If at least one of the indices is 1,
% %                       % then it is an open region and we can't
% %                       % patch that.
% %         %polygon = [polygon; V1(C{i})];
% %         %patch(V1(C{i},1),V1(C{i},2),i); % use color i.
% %     end
% %     if ~isempty(C2{i})
% %         %patch(V2(C2{i},1),V2(C2{i},2),i);
% %     end
%     patch(polygon(C_final{i},1),polygon(C_final{i},2),i);
% end




%% JUNK/TEMPORARY CODE:

%[B,I] = sort(minInd);

% Cbb = [bbii(:,1) bbii(:,1)+1]
% %bbox(bbii(:,1),:)
% [Cbb bbx bby] % V2 extra
% 
% test = Cbb(:,1)
% 
% C3 = num2cell(repmat(NaN,1,length(Cbb))');
% for i=1:length(test)
%     C3{i} = find(minInd==test(i));
% end
%     
% C3 = num2cell(repmat(NaN,1,length(Cbb))');
% for i=1:length(C2)
%    C2{i}
%    %C3{i} = minInd(Cbb(i,:));
% end
% 
% minInd


% bbx2 = [bbx; bbox(1:end-1,1)];
% bby2 = [bby; bbox(1:end-1,2)];


% % loop though all points with open regions and select boundary points
% n = length(bbx2);
% for i = 1:length(C)
%     if any(C{i}==1) % If at least one of the indices is 1,
%         % then it is an open region and we can't patch that.
%         linex = [bbx2'; repmat(xy(i,1), 1,n)];
%         liney = [bby2'; repmat(xy(i,2), 1,n)];
%         [xi_temp, yi_temp, ii] = polyxpoly(linex, liney, vx, vy);
% %         for j = 1:length(linex)
% %             [xi, yi] = polyxpoly(vx(:,j), vy(:,j), bbox(:,1), bbox(:,2));
% %             bbx = [bbx; xi];
% %             bby = [bby; yi];
% %         end
%     end
% end

figure(2);
plot(vx_extent,vy_extent);
xlim([-6 6]); ylim([-6 6]);
hold on;
plot(bbox(:,1),bbox(:,2));
scatter(bbx, bby);
scatter(V1(:,1), V1(:,2));
%scatter(bbox(1:end-1,1),bbox(1:end-1,2));
text(bbox(1:end-1,1)-0.1, bbox(1:end-1,2)-0.1, cellstr(num2str([1:4]')) );
scatter( xy(:,1)' , xy(:,2)', 'filled');
text(xy(:,1)'+0.1, xy(:,2)'+0.1, cellstr(num2str([1:n]')) );
%plot(linex,liney,':');
%scatter(xi_temp, yi_temp);
%patch(test(:,1),test(:,2),3);
hold off;
% xlabel('x'); ylabel('y');
% %voronoi(xy(:,1)' , xy(:,2)');
% hold off;

%[xi, yi] = polyxpoly(vx, vy, bbox(:,1), bbox(:,2),'unique'); % line intersect with boundary box

%[lat,lon] = polymerge(vx,vy);
%plot(lat,lon, ':');

%[XY,V1,C] = VoronoiLimit(xy(:,1)' , xy(:,2)', bbox); % http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit


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

figure(5);
plot(bvx_split,bvy_split);
xlim([-6 6]); ylim([-6 6]);

figure(33);
plot(vx_extent,vy_extent);
xlim([-8 8]); ylim([-8 8]);
hold on;
plot(bbox(:,1),bbox(:,2));
scatter(out.intMatrixX(out.intAdjacencyMatrix),out.intMatrixY(out.intAdjacencyMatrix),[],'r');
hold off;



% Calculating distance from line segments to points

bdist = pldist2(xy, BXY)
% Find closest firm for each line segment
[~, I] = min(bdist)



