function D = pldist2(points, lines)
%PLDIST2 calculates the 2D point-line distance between every point and every
% line. The input POINTS is a matrix where every row is a point with form
% [x y]. And the inpit LINES is a matrix where every row is a line with
% form [x1 y1 x2 y2].
%
% Example:
%    points = [0 0;
%              1 1]; % x y
%    lines = [1 2 3 5;
%             2 2 5 7]; % x1 y1 x2 y2
%   
%   Jonas K. Sekamane. 
%   Version 0.04
%
%   Inspired in part by: 
%       http://stackoverflow.com/a/25803408/1053612

dist = NaN( size(points,1), size(lines,1) );
for l=1:size(lines,1)
    a = lines(l,1:2); % segment points a,b
    b = lines(l,3:4);
    
    for p=1:size(points,1)
        x = points(p,:);
        
        d_ab = norm(a-b);
        d_ax = norm(a-x);
        d_bx = norm(b-x);

        if dot(a-b, x-b) * dot(b-a, x-a) >= 0
            A = [a 1; b 1; x 1];
            dist(p,l) = abs(det(A)) / d_ab;     
        else
            dist(p,l) = min(d_ax, d_bx);
        end
    end
end

D = dist;

end