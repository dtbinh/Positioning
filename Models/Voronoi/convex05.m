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

% Triangle, Point, Angle
t_points = DT.ConnectivityList';
t_angles = tri_angle';
[repelem(1:n_tri, 3)' t_points(:) t_angles(:)]


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
