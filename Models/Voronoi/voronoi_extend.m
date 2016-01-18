function [vx_extend, vy_extend] = voronoi_extend(vx, vy, bbox, bdistance)
%VORONOI_EXTEND
%   Extends the lines to insure that they intersect with the boundary box.
%   Only extend lines where the endspoints are unconnected with other lines 
%   (and extend in the direction away from the connected end). The length 
%   after extension should be the diagonal length of the boundary box.
%   
%   Jonas K. Sekamane. 
%   Version 0.01
%

	% Identify the Voronoi edges/endpoints that are within the boundary.
	[~, vui] = unique(vx(:)); % Only interested in unique endpoints.
	V0 = [vx(sort(vui)) vy(sort(vui))]; % Set of unique xy-endpoints.
	in = inpolygon(V0(:,1), V0(:,2), bbox(:,1), bbox(:,2)); % Index endpoins within boundary.


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
    vx_ext = vx;
    vy_ext = vy;
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
        vx_ext(j,line) = vx(i,line) + dx;
        vy_ext(j,line) = vy(i,line) + dy;
    end

    % Output variables
    vx_extend = vx_ext;
    vy_extend = vy_ext;
end