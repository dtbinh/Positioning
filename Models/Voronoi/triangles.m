function [convexTriangleSets,DT] = triangles(xy)
    % inputs
    %xy = [3.3735 0.7889; -0.1072 -3.4814; -3.9732 4.1955; -5 5; 5 5; 5 -5; -5 -5];
    DT = delaunayTriangulation(xy);

    function convexTriangleSets = testAddTriangle(iTriangle,attachedTriangleIndices,includedTriangleIndices)
        % add triangle
        includedTriangleIndices(end+1) = iTriangle;

        % naive test if allpoint of set of triangles are on convex hull
        nodes = unique(DT(includedTriangleIndices,:));
        coords = DT.Points(nodes,:);
        ch = convexHull(delaunayTriangulation(coords)); 
        allNodesOnConvexHull = length(nodes) == length(ch)-1;

        if ~allNodesOnConvexHull
            convexTriangleSets = {};
            return 
        end

        % find triangles connected to iTriangle
        currentTriangle = DT.ConnectivityList(iTriangle,:)';

        attachedCell = DT.edgeAttachments(currentTriangle([1 2 3]),currentTriangle([2 3 1]));
        attachedRow = unique([attachedTriangleIndices,attachedCell{:}]);
        attachedTriangleIndices = attachedRow(~ismember(attachedRow,includedTriangleIndices));

        % recursively try to expand connected triangles
        convexTriangleSets = {sort(includedTriangleIndices)};
        for ii = 1:length(attachedTriangleIndices)
            convexTriangleSets = [convexTriangleSets,...
                testAddTriangle(attachedTriangleIndices(ii),...
                                attachedTriangleIndices,...
                                includedTriangleIndices)]; %#ok<AGROW>
        end
    end

    includedTriangleIndices = [];
    attachedTriangleIndices = [];
    convexTriangleSets = {};
    for iTriangle = 1:DT.size
        convexTriangleSets = [convexTriangleSets,...
            testAddTriangle(iTriangle,attachedTriangleIndices,includedTriangleIndices)]; %#ok<AGROW>
    end

    % filter single triangles
    %convexTriangleSets(cellfun(@length,convexTriangleSets) == 1) = [];

    % filter unique sets; convert to string because matlab cannot unique a cell array
    [~,c] = unique(cellfun(@(x) sprintf('%d,',x),convexTriangleSets,'UniformOutput',false));
    convexTriangleSets = convexTriangleSets(c);
    
    % Get all the convex sets with at least two triangles.
    [L, I] = sort(cellfun(@length,convexTriangleSets),'descend');
    L_idx = I(L>1);

    combi = cell(0);
    % Loop through these convex sets.
    for ll = 1:length(L_idx)
        % Extract array with triangles
        A = convexTriangleSets{L_idx(ll)};
        % Calculate all using these triangles. (don't include full array itself).
        for i = 1:size(A,2)-1
            combi = [combi; num2cell(nchoosek(A,i),2)];
        end
    end
    % Remove duplicates
    [~,c] = unique(cellfun(@(x) sprintf('%d,',x),combi,'UniformOutput',false));
    combi = combi(c);

    % Remove the convex sets that are not the largest.
    rm = ismember(cellfun(@(x) sprintf('%d,',x),convexTriangleSets','UniformOutput',false), cellfun(@(x) sprintf('%d,',x),combi,'UniformOutput',false));
    convexTriangleSets(rm) = [];
    

    % plot result
    n = ceil(sqrt(length(convexTriangleSets)));
    for kk = 1:length(convexTriangleSets)
        subplot(n,n,kk)
        triplot(DT,'k')
        hold on
        patch('faces',DT(convexTriangleSets{kk},:), 'vertices', DT.Points, 'FaceColor','r');
    end
end