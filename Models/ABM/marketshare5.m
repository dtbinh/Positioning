function [M, U] = marketshare5(firms, customers)
%MARKETSHARE5   calculate the marketshare of each firm
%   MARKETSHARE5(FIRMS, CUSTOMERS) 
%   where CUSTOMERS is a matrix of all the coordinates of the market and
%   FIRMS contains the coordinates of the firm.
%
%   Output variables:
%   * M is a matrix of the entire market where the value of each cell is the closest firm.
%   * U is a matrix of the disutility in the entire market.
%   
%   Jonas K. Sekamane. 
%   Version 0.05

    %% Setup
    
    %% Closest firm
    
    if(size(firms,1)>1)
        [all, C, hull] = voronoib(firms);
        
        % Determine market based on convex hull / polygons
        market = NaN( length(customers), 1);
        for i = 1:size(firms,1)
             X = all(C{i},1);
             Y = all(C{i},2);
             idx = inpolygon(customers(:,1), customers(:,2), X(hull{i}), Y(hull{i}));
             market(idx) = i;
             % TO-DO: inpolygon ON the boundary. Then randomly pick firm
        end
        
        minInd = market; 
        
    else
        minInd = ones(length(customers), 1);
    end

    %% Shape output variables
    l = sqrt(length(customers));
    customers = reshape(minInd,[l l]);
    
    % Output variables
    M = customers;
    U = NaN(l,l);
end