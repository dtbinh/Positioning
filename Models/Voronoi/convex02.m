xy = [3.3735 0.7889; -0.1072 -3.4814; -3.9732 4.1955; -5 5; 5 5; 5 -5; -5 -5];
DT = delaunayTriangulation(xy);
%The coordinates for each triangle -- each row is a triangle.
TRIX = reshape(DT.Points(DT.ConnectivityList, 1), size(DT.ConnectivityList));
TRIY = reshape(DT.Points(DT.ConnectivityList, 2), size(DT.ConnectivityList));
figure(20);
triplot(DT, 'k');
hold on;
scatter(xy(:,1), xy(:,2), 'filled', 'b');
text(xy(:,1), xy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','Blue','FontSize',14);
text(mean(TRIX, 2), mean(TRIY, 2), cellstr(num2str([1:size(TRIY,1)]')), 'Color','Red');
xlim([-6 6]); ylim([-6 6]);
set(gca,'YTick',(-5:5:5));
set(gca,'XTick',(-5:5:5));


% Difference between points
diffx = diff([TRIX TRIX(:,1)], [], 2);
diffy = diff([TRIY TRIY(:,1)], [], 2);
diffxy = [diffx(:) diffy(:)];
% Norm
normxy = reshape( arrayfun(@(row) norm(diffxy(row,:)), 1:size(diffxy,1)), size(DT.ConnectivityList));
% sik^2 + sjk^2 - sij^2, for i,j,k = 1,2,3 and i?j?k.
nominator = repmat(sum(normxy.^2, 2), 1, 3) - 2*normxy.^2;
% 2 * sik * sjk , for i,j,k = 1,2,3 and i?j?k.
denominator = 2 * repmat(prod(normxy, 2), 1, 3)./normxy;
% angle
tri_angle = acosd(nominator./denominator);
tri_angle = circshift(tri_angle, [0 -1]);

n_tri = size(TRIX,1);

% Adjacency matrix connecting points (rows) with triangles (columns).
adj_points = zeros(size(xy,1), n_tri);
adj_angle = NaN(size(adj_points));
for point =1:size(xy,1)
    idx = find(DT.ConnectivityList == point);
    [a_tri, ~] = ind2sub(size(DT.ConnectivityList), idx);
    adj_points(point,a_tri) = 1;
    adj_angle(point,a_tri) = tri_angle(idx);
end
adj_points = logical(adj_points);

% All edges in the Delaunay triangulation
DT_edges = edges(DT);
% Adjacency connecting edges (rows) with triangles (columns).
adj_edge = logical(adj_points(DT_edges(:,1),:) .* adj_points(DT_edges(:,2),:));

edgesangles = NaN(size(DT_edges));
adj = zeros(n_tri);
adj_convex = zeros(n_tri);
for edge=1:size(DT_edges,1)
    % The angles on either side of the edge.
    tri = adj_edge(edge,:);
    t = adj_angle(DT_edges(edge,:), tri );
    edgesangles(edge,:) = sum(t, 2);
    tri_idx = find(tri);
    adj(tri_idx,tri_idx) = 1;
    adj_convex(tri_idx,tri_idx) = prod(edgesangles(edge,:) <= 180);
end
convexedges = (edgesangles <= 180);
% Set diagonals to zero.
adj(logical(eye(n_tri))) = 0;
adj_convex(logical(eye(n_tri))) = 0;

%ec = logical(prod(convexedges,2));
%[DT_edges(ec,:) edgesangles(ec,:)]







%%
adj_convex2 = adj_convex^2;
adj_convex2(logical(eye(n_tri))) = 0;

[tri1, tri2] = find(triu(adj_convex2));

% The triangles common point and common triangle
[~, tri_points] = max( adj_points(:,tri1) .* adj_points(:,tri2) , [], 1);
[~, tri_common] = max( adj(tri1,:) .* adj(tri2,:), [], 2)
[~, tri_edge1] = max( adj_edge(:,tri1) .* adj_edge(:,tri_common) , [], 1);
[~, tri_edge2] = max( adj_edge(:,tri2) .* adj_edge(:,tri_common) , [], 1);

index = ([DT_edges(tri_edge1, :) DT_edges(tri_edge2, :) ]-repmat(tri_points', 1, 4) == 0)
angles_common = [edgesangles(tri_edge1, :) edgesangles(tri_edge2, :)]'
angles_common_point = reshape(angles_common(index'), 2, size(tri_points,2))';

tri_convex = (sum(angles_common_point,2) <= 180);

[tri1 tri2 tri_convex]


adj_convex_temp = adj_convex;
rm_idx = sub2ind(size(adj_convex), [tri1(tri_convex); tri2(tri_convex); repmat(tri_common(tri_convex),2,1)], [repmat(tri_common(tri_convex),2,1); tri1(tri_convex); tri2(tri_convex)] )
adj_convex_temp(rm_idx) = 0;
[tria, trib] = find(triu(adj_convex_temp));



F1 = find(sum(adj_convex, 2)==0);
F2 = [tria trib];
F3 = [tri1(tri_convex) tri2(tri_convex) tri_common(tri_convex)];

cs = cumsum([size(F1,1) size(F2,1) size(F3,1)]);

convexhull_matrix = NaN( cs(end), 3 );
convexhull_matrix(1:cs(1), 1) = F1;
convexhull_matrix(cs(1)+1:cs(2), 1:2) = F2;
convexhull_matrix(cs(2)+1:cs(3), 1:3) = F3;

convexhull = num2cell(convexhull_matrix, 2);
%convexhull(cellfun(@isnan, convexhull, 'UniformOutput', false)) = [];


%[tri1 tri2 tri_common tri_points' tri_edge1' tri_edge2']

%FINAL_adj_convex = adj_convex;
%FINAL_adj_convex(tri1(tri_convex),tri2(tri_convex)) = 1;
%triu(FINAL_adj_convex)




adj_edge(:,tri1)  .* adj_edge(:,tri2) .* adj_edge(:,tri_common)



for i = 1:size(tri1,1)
    adj(tri1(i),:) .* adj(tri2(i),:)
    
end    
%adj_points(sub2ind(size(adj_points), tri_points', tri1))


adj_angle(:,tri1) + adj_angle(:,tri2)

[Loca, Locb] = ismember([tri1 tri2], DT_edges, 'rows')
% 2, 6, NaN, NaN, 10

DT_edges(Locb,:)

[DT_edges edgesangles]

for point = 1:size(xy,1)
    point_edge = vertexAttachments(DT, point)
    
    edgesangles(point_edge{1},:)
    DT_edges(point_edge{1},:)
    
    %A_tri = find(adj_points(point,:));
    %A = adj_angle(point,A_tri);
    %C = cell(0);
    %for i = 1:size(A,2)
    %    C = [C; num2cell(nchoosek(A,i),2)];
    %end
    %C = C(cellfun(@(x) sum(x), C) <= 180);
end
