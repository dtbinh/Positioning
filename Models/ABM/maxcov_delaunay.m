function xy_new = maxcov_delaunay( firms, customers, F )
%MAXCOV_DELAUNAY
%   Determines the position of a new firm entering a market with exisiting
%   firms, that maximises its market share.
%
%   Returns the coordinates of the weighted centroid of the triangle with 
%   the largest market share. The triangles are constructed using the 
%   delaunay triangulation from the position of existing firms and the 
%   boundary points. The method takes into account the underlying 
%   population density.
%   
%   Jonas K. Sekamane. 
%   Version 0.01
   
    % Boundary box
    bbox = [-5 5 5 -5 -5; 5 5 -5 -5 5]';

    % Adding the coordinates of the boundary points to the list of firm
    % coordinates. Allowing for the new firm's position to be outside the
    % convex hull created by exsisting firm positions.
    xy_boundary = [unique(firms, 'rows'); bbox(1:end-1,1) bbox(1:end-1,2)];
    % Create Delaunay Triangulation
    DT = delaunayTriangulation(xy_boundary);
    % Count the number of triangles in the Delaunay Triangulation.
    triangles = size(DT.ConnectivityList, 1);

    % The XY coordinates of the triangle
    DT2X = reshape(DT.Points(DT.ConnectivityList, 1), size(DT.ConnectivityList));
    DT2Y = reshape(DT.Points(DT.ConnectivityList, 2), size(DT.ConnectivityList));
    
    F_tri  = NaN(triangles,1);
    centroid_tri  = NaN(triangles,2);
    % Loop through all triangles
    for tri = 1:triangles
        
        % Index of customers within the triangle
        idx = InPolygon(customers(:,1), customers(:,2), DT2X(tri,:), DT2Y(tri,:));
        
        % Number of customers in the triabgle
        F_tri(tri) = sum( F(idx) );
        
        % The xy-coordinate of the centroid within the triangle weighted with the probability density
        centroid_tri(tri,:) = ([customers(idx,1) customers(idx,2)]' * F(idx))' ./ F_tri(tri);
        
    end
    
    % Market share of each triangle.
    %F_tri/sum(F(:))
    
    % Find the index of the triangle with the largest market share.
    [~, max_tri] = max(F_tri/sum(F(:)));

    % Centroid oordinates of the triangle with largest market share.
    %centroid_tri(max_tri,:)
    
    
    % Output variables
    xy_new = centroid_tri(max_tri,:);
end