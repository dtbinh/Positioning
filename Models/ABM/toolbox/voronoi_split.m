function [VXY_split, BXY_split] = voronoi_split(VXY, BXY, A, IX, IY)
%VORONOI_SPLIT
%   Splits line segments into new line segments at the intersection points. 
%   The two sets of lines are given by VXY and BXY. In both sets each row 
%   is a line with form [x1 y1 x2 y2].
%   The adjacency matrix A indicates which lines intersect. The columns 
%   correspond to the lines in my boundary set (BXY), and the rows 
%   correspond to the lines in my Voronoi set (VXY). The matricies IX and
%   IY contain the coordinates for the intersection points. The two
%   matricies are formated similar to the adjacency matrix.
%   
%   Jonas K. Sekamane. 
%   Version 0.01
%

    % New variables that simply later code
    %A = out.intAdjacencyMatrix;
    %IX = out.intMatrixX;
    %IY = out.intMatrixY;
    s = size(A);

    % Create new matrix for the split line segments
    VXY_s = NaN(s(1)+sum(A(:)), 4);
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
        VXY_s(v:v-1+l,:) = reshape(pointssorted', 4, l)';
        v = v+l;

    end
    VXY_s(v:end,:) = VXY(vbreaks==0,:); % Lines with no split

    % Create new matrix for the split line segments
    BXY_s = NaN(s(2)+sum(A(:)), 4);
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
        BXY_s(b:b-1+l,:) = reshape(pointssorted', 4, l)';
        b = b+l;

    end
    BXY_s(b:end,:) = BXY(bbreaks==0,:); % Lines with no split

    % Delete lines outside boundary (.0000000000001 fixes rounding error).
    %VXY_s(find( sum( (abs(VXY_s) > abs(BXY(1,1))+.0000000000001) ,2) ),:) = [];
    
    
    % Output variables
    VXY_split = VXY_s;
    BXY_split = BXY_s;
end