function F = marketsharen(firms, customers)
%MARKETSHAREN   calculate the marketshare of each firm
%   MARKETSHAREN(FIRMS, CUSTOMERS) 
%   where CUSTOMERS is a graph/adjacency matrix of all the customers and
%   and FIRMS contains the number of the firm in the graph lattice.
%   Jonas K. Sekamane. 
%   Version 0.01

    % Array with each customer and their closest firm.
    for n = 1:length(firms)
       [dist, path, pred] = graphshortestpath(customers, firms(n), 'Directed', false);
       distance(n,:) = dist;
       %ShortPath(n,:) = path;
       %ShortPred(n,:) = pred;
    end

    % Find for each customer the minimum distiance (minVal) and the firm (minInd)
    [minVal, minInd] = min(distance);

    % A customers with equal distiance to two firms will randomly pick one firm.
    % Loop through all distances and find rows that equal the minimum value.
    % If there is more than 1 row, then randomly draw.
    for i = 1:size(distance,2)
        rows = find(distance(:,i)==minVal(i));
        if size(rows,1) > 1
            minInd(i) = randsample(rows,1);
        end
    end
    
    % Array with each customer/node and their closest firm.
    F = minInd;
end