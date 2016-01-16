%% Voronoi diagram
% version 0.01
% Jonas K. Sekamane
%
% *BOUNDED VORONOI DIAGRAM*
% Extending the the _voronoi_ and _voronoin_ functions with a boundary, 
% such that any open regions can be patched.

clearvars;

% Points - Random firm locations
n = 4;
x0 = rand(1, n)*10-5; % Initial x-position of firm
y0 = rand(1, n)*10-5; % Initial y-position of firm
xy = [x0' y0'];

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


% Identify Voronoi line segments with common endpoints.
vend = NaN(size(vx));
% for i=1:2
%     vend(i,:) = ismember(vx(i,:), V0(:,1));
% end
for line=1:length(vx);
    for i=1:2
        vend(i,line) = ismember(vx(i,line), vx(:,1:end~=line));
    end
end
% Find the Voronoi lines that have exactly one shared endpoint
count = sum(vend,1);
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


% Line intersect with boundary box
bbx = []; bby = []; bbii = [];
for j = 1:length(vx)
    [xi, yi ii] = polyxpoly(bbox(:,1), bbox(:,2), vx(:,j), vy(:,j)); % [xi, yi, ii] = 
    bbx = [bbx; xi];
    bby = [bby; yi];
    bbii = [bbii; ii];
end
V3 = [bbx bby]; % Coordinates of the voronoi diagram lines intersection with the boundary box.

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
plot(vx,vy);
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


figure(3);
plot(vx,vy);
xlim([-30 30]); ylim([-30 30]);
hold on;
plot(bbox(:,1),bbox(:,2));
hold off;

figure(31);
plot(vx_extent,vy_extent);
xlim([-30 30]); ylim([-30 30]);
hold on;
plot(bbox(:,1),bbox(:,2));
hold off;

% figure(4);
% plot(vx(:,bbii(1,2)),vy(:,bbii(1,2)));
% xlim([-6 6]); ylim([-6 6]);
% hold on;
% plot( bbox(bbii(1,1),1) , bbox(bbii(1,1),2) );
% hold off;
% bbii(1,1)
% 
% figure(5);
% plot(bvx,bvy);
% xlim([-6 6]); ylim([-6 6]);

