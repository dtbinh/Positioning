function F = stairs2(y, varargin)
%STAIRS2 is similar to STAIS, but only displays the vertical lines.
%   Jonas K. Sekamane. 
%   Version 0.01
%   Inspired by http://stackoverflow.com/a/27337323/1053612
   
        for i=1:length(y)-1
            A(:,i) = [i (i+1)];
            B(:,i) = [y(i) y(i)];
        end
        plot(A,B,varargin{:});
end