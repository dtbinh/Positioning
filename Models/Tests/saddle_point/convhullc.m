function C = convhullc(pp, ll)
% Find convex hulls around a set of points, where the convex hulls are
% constrained by a set line segments.
% In the pp matrix each row is a point [x y].
% In the ll matrix each row is a line with format [x1 y1 x2 y2].
    
    np = size(pp,1); % Number of points
    nl = size(ll,1); % Number of lines
    
    % create line connecting all pairs of points
    combi = nchoosek(1:np, 2);
    ppxy = [pp(combi(:,1),:) pp(combi(:,2),:)];
    
    % Intersection between all lines connecting pairs of points and all line segments
    intersectpl = lineSegmentIntersect( ppxy, ll );
    % Pairs of points that do not intersect a line segment
    nointersect = (sum(intersectpl.intAdjacencyMatrix,2) == 0);
    % combi(nointersect,:)
    % Creating adjacency matrix of points with no intersect 
    adj = sparse(combi(nointersect,1), combi(nointersect,2), 1, np, np);
    % Create undirected graph.
    adj = adj + adj.';
    % Find the connected components
    [nch, bin] = conncomp(adj);

    ch = cell(nch,1);
    for i=1:nch
        if( sum((bin==i)) > 2 )
            ppx = pp((bin==i),1);
            ppy = pp((bin==i),2);
            if (size(ppx,1) > 3)
                K = convhull(ppx, ppy);
                ch{i} = [ppx(K) ppy(K)];
            else
                ch{i} = [];
            end
        end
    end
    ch = ch(~cellfun('isempty',ch));

    C = ch;
end