function F = regionpropsext(varargin)
%REGIONPROPSEXT Extends REGIONPROPS. 
%   Measures the perimeter of the extrema points.
%   REGIONPROPSEXT(BW,PROPERTIES) 
%   where PROPERTIES must include 'Extrema'.
%   Jonas K. Sekamane. 
%   Version 0.01

    share_props = regionprops(varargin{:});
    
    for n = 1:size(share_props,1)
        share_props(n).Extrema = [share_props(n).Extrema; share_props(n).Extrema(1,:)]; % Close the loop/polygon by "connecting" last and first extrema point.
        % Square the differences in respectively the x-coordinates and the y-coordinates;
        % Sum the squared terms and take the square root to get the distance between each extrema point;
        % Finally sum over all distances to get total distance or perimeter;
        share_props(n).ExtremaPerimeter = sum( sqrt(sum( diff(share_props(n).Extrema).^2 , 2)) ); % Euclidean distance
    end
    
    F = share_props;