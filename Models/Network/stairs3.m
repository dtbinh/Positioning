function F = stairs2(y, varargin)
%STAIRS2 is similar to STAIS, but horizontal lines are dimmed.
%   STAIRS2(y, 'Alpha', 0.1, 'Color', [0 0 0], ...)
%   where ALPHA input must be first followed by COLOR input.
%   Jonas K. Sekamane. 
%   Version 0.01
%   Inspired by http://stackoverflow.com/a/27337323/1053612
        
        alpha = varargin{2};
        color = varargin{4};

        for i=1:length(y)-1
            A(:,i) = [i (i+1)];
            B(:,i) = [y(i) y(i)];
            
            C(:,i) = [i+1 i+1];
            D(:,i) = [y(i) y(i+1)];
        end
        plot(A,B, 'Color', color, varargin{5:end});
        hold on;
        plot(C,D, 'Color', [color alpha], varargin{5:end});
end