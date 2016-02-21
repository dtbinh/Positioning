
xy0 = [3.3735 0.7889; ...
      -0.1072 -3.4814; ...
      -3.9732 4.1955; ...
      -5 5; ...
       5 5; ...
       5 -5; ...
      -5 -5];
xy1 = [-5 5; ...
       5 5; ...
       3.3735 0.7889; ...
       -3.9732 4.1955];
xy2 = [5 5; ...
       5 -5; ...
       3.3735 0.7889];
xy3 = [3.3735 0.7889; ...
      -3.9732 4.1955; ...
      -0.1072 -3.4814; ...
       5 -5];
xy4 = [-0.1072 -3.4814; ...
       5 -5; ...
      -5 -5];
xy5 = [-0.1072 -3.4814; ...
      -3.9732 4.1955; ...
      -5 5; 
      -5 -5];
xy6 = [-0.1072 -3.4814; ...
       3.3735 0.7889; ...
       5 5; ...
      -3.9732 4.1955; ... 
      -5 -5];
xy7 = [3.3735 0.7889; ...
      -0.1072 -3.4814; ...
      -3.9732 4.1955; ...
      -5 5; 
       5 5; ...
       5 -5; ...
      -5 -5];
xy8 = [3.3735 0.7889; ...
      -0.1072 -3.4814; ...
      -3.9732 4.1955; ...
      -5 5; 
       5 5; ...
       5 -5; ...
      -5 -5];


 
 
 figure(12);
 clf reset; % Reset figure.
 scatter(xy0(:,1), xy0(:,2));
 xlim([-6 6]); ylim([-6 6]);
 set(gca,'YTick',(-5:5:5));
 set(gca,'XTick',(-5:5:5));
 hold on;
 %patch(xy1(:,1), xy1(:,2), 1, 'FaceAlpha', 0.2);
 %patch(xy2(:,1), xy2(:,2), 2, 'FaceAlpha', 0.2);
 %patch(xy3(:,1), xy3(:,2), 3, 'FaceAlpha', 0.2);
 %patch(xy4(:,1), xy4(:,2), 4, 'FaceAlpha', 0.2);
 %patch(xy5(:,1), xy5(:,2), 5, 'FaceAlpha', 0.2);
 patch(xy6(:,1), xy6(:,2), 6, 'FaceAlpha', 0.2);
 hold off;
 
 
 %% 
 
xy = [3.3735 0.7889; -0.1072 -3.4814; -3.9732 4.1955; -5 5; 5 5; 5 -5; -5 -5];
DT = delaunayTriangulation(xy);

figure(20);
triplot(DT);
text(xy(:,1), xy(:,2)+0.25, cellstr(num2str([1:size(xy,1)]')), 'Color','red','FontSize',14);
text(mean(TRIX, 2), mean(TRIY, 2), cellstr(num2str([1:size(TRIY,1)]')), 'Color','blue');



%The coordinates for each triangle -- each row is a triangle.
TRIX = reshape(DT.Points(DT.ConnectivityList, 1), size(DT.ConnectivityList));
TRIY = reshape(DT.Points(DT.ConnectivityList, 2), size(DT.ConnectivityList));

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
angle = acosd(nominator./denominator);
angle = circshift(angle, [0 -1]);

n_tri = size(TRIX,1);



% Adjacency connecting points (rows) with triangles (columns).
adjacency_points = zeros(size(xy,1), n_tri);
adjacency_angle = NaN(size(adjacency_points));
for point =1:size(xy,1)
    idx = find(DT.ConnectivityList == point);
    [a_tri, ~] = ind2sub(size(DT.ConnectivityList), idx);
    adjacency_points(point,a_tri) = 1;
    adjacency_angle(point,a_tri) = angle(idx);
end


% All edges in the Delaunay triangulation
edgespoints = edges(DT);
% Adjacency connecting edges (rows) with triangles (columns).
adjacency_edge = logical(adjacency_points(edgespoints(:,1),:) .* adjacency_points(edgespoints(:,2),:));

edgesangles = NaN(size(edgespoints));
adj = zeros(n_tri);
adjacency_convex = zeros(n_tri);
for edge=1:size(edgespoints,1)
    % The angles on either side of the edge.
    tri = adjacency_edge(edge,:);
    t = adjacency_angle(edgespoints(edge,:), tri );
    edgesangles(edge,:) = sum(t, 2);
    tri_idx = find(tri);
    adj(tri_idx,tri_idx) = 1;
    adjacency_convex(tri_idx,tri_idx) = prod(edgesangles(edge,:) <= 180);
end
convexedges = (edgesangles <= 180);
% Set diagonals to zero.
adj(logical(eye(n_tri))) = 0;
adjacency_convex(logical(eye(n_tri))) = 0;


active_e = edgespoints;
while ~isempty(active_e)
    for edge=1:size(active_e,1)
        tri = adjacency_edge(edge,:);
        t = adjacency_angle(active_e(edge,:), tri );
        prod( (sum(t, 2) <= 180) )
    end
    
    active_e = [];
end




prod(convexedges,2)

[edgespoints edgesangles]

adjacency_convex

adjacency_convex2 = adjacency_convex^2;
adjacency_convex2(logical(eye(n_tri))) = 0;

adjacency_convex3 = adjacency_convex2*adjacency_convex;
adjacency_convex3(logical(eye(n_tri))) = 0;
adjacency_convex3*adjacency_convex

prod(convexedges,2)

(adj^2 > 0)

(adj^5 > 0)

size(adjacency_edge')

sum(adjacency_edge, 2)


% t(1, [1 2]) + t(2, [2 1])
% 
% reshape(t(:), [], 2)
% flipud(t)
% trace(t)
% rot90(t)
% angle_edge(edge,:)
    
%sum(test, 2, 'omitnan')

%[a b]=find(DT.ConnectivityList==1);
%[c d]=find(DT.ConnectivityList==2);
%intersect(a,c)

%[row_edge, col_tri] = find( adjacency(edgespoints(:,1),:) .* adjacency(edgespoints(:,2),:) )



[a_edge, a_tri] = find( adjacency_edge );
[~, a_order] = sort(a_edge);
EdgesList = [a_edge(a_order) a_tri(a_order)];





adjacency_angle(edgespoints(:,1), EdgesList)

test1 = adjacency_angle(edgespoints(:,1),:);
test2 = adjacency_angle(edgespoints(:,2),:);


adjacency_edge2 = (adjacency_edge .* cumsum(adjacency_edge, 2) == 2);
adjacency_edge1 = adjacency_edge - adjacency_edge2;

sumangle = sum(test1 .* adjacency_edge1 + test2 .* adjacency_edge2, 2, 'omitnan')
edgespoints((sumangle > 180),:)



[m where] = max( cumsum(adjacency_edge, 2) ,[], 2);
where = where .* (m==2)

adjacency_edge2 = zeros(size(adjacency_edge));

adjacency_edge2( sub2ind(size(adjacency_edge), where, 1:size(adjacency_edge,1)) ) = 1

arrayfun(@(i) find(adjacency_edge2(i,:)==2,1,'first'),1:size(adjacency_edge2,2))

adjacency_edge2 = (cumsum(adjacency_edge, 2) == 2);
[dum,idx] = sort(adjacency_edge2,2,'descend');
firstOne = idx(:,1)
adjacency_edge = (cumsum(adjacency_edge, 2) == 2)

test1(adjacency_edge)
test2(adjacency_edge)


%test = test1 + test2;


size(col)
reshape(col, [], 2)


adjacency_points(edgespoints(:,1),:) .* adjacency_points(edgespoints(:,2),:)


ismember([1 2], DT.ConnectivityList)
