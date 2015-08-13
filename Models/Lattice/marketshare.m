function F = marketshare(firms, customers)
%MARKETSHARE   calculate the marketshare of each firm
%   MARKETSHARE(FIRMS, CUSTOMERS) 
%   where CUSTOMERS is a lattice of all the customers and
%   and FIRMS contains the coordinates of the firm in this lattice.
%   Jonas K. Sekamane. 
%   Version 0.01

    % Create empty matrix for coordinates. First coloumn is x and second is y
    coordinates = cell(numel(customers),2);
    for idx = 1:numel(customers)
        I = cell(1, ndims(customers));
        [I{:}] = ind2sub(size(customers),idx);
        coordinates(idx,:) = I(:,[2,1]);
    end
    coordinates = cell2mat(coordinates);

    % For each customer calculate the euclidean distance to each of the firms.
    distance = pdist2(firms,coordinates,'euclidean');

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
    
    % Fills out the customer lattice with the ID of the closest firm.
    for idx = 1:numel(customers)
        customers(idx) = minInd(idx);
    end;
    
    % Matrix with each customer and their closest firm.
    F = customers;
end