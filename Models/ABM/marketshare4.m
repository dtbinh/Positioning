function [M, U] = marketshare4(firms, customers)
%MARKETSHARE4   calculate the marketshare of each firm
%   MARKETSHARE4(FIRMS, CUSTOMERS) 
%   where CUSTOMERS is a matrix of all the coordinates of the market and
%   FIRMS contains the coordinates of the firm.
%
%   Output variables:
%   * M is a matrix of the entire market where the value of each cell is the closest firm.
%   * U is a matrix of the disutility in the entire market.
%   
%   Jonas K. Sekamane. 
%   Version 0.04

    %% Setup

    % Create empty matrix for coordinates. First coloumn is x and second is y
%     coordinates = cell(numel(customers),2);
%     for idx = 1:numel(customers)
%         I = cell(1, ndims(customers));
%         [I{:}] = ind2sub(size(customers),idx);
%         coordinates(idx,:) = I(:,[2,1]);
%     end
%     coordinates = cell2mat(coordinates);
    
    %% Closest firm
    
    % For each customer calculate the euclidean distance to each of the firms.
    distance = pdist2(firms, customers, 'euclidean');

    if(size(firms,1)>1)
    
        % Find for each customer the minimum distiance (minVal) and the firm (minInd)
        [minVal, minInd] = min(distance);
        
        % The absolute difference between the distances and the minimum
        % values.
        absdiff = abs(distance-repmat(minVal, length(firms),1));
        % If the difference is zero then this is minimum distance.
        minValIdx = absdiff==0; 
        % Check if customers has equal distiance to two or more firms.
        c = find(sum(minValIdx,1)>1);
        % Loop through all customers with equal distiance to several firms 
        % and randomly pick one firm.
        for idx = c
            rows = find(minValIdx(:,idx)==1);
            minInd(idx) = rows(randi(length(rows))); % use rows(randi(length(rows))) rather than randsample(rows,1) since it relies on seed from rnd()
        end
         
    else
        minInd = ones(size(distance));
        minVal = distance;
    end
    
    % Fills out the customer lattice with the ID of the closest firm.
%     for idx = 1:numel(customers)
%         customers(idx) = minInd(idx);
%     end;
         
    %% Centroid
    
    % The centroid coordinates for each market
%     centroid = cell2mat( ...
%                     arrayfun(@(firm) ...
%                         mean(customers(find(minInd==firm),:),1), ... % Calculate mean xy-coordinates of each firms market
%                         1:length(firms), ... % Firm ID/name
%                         'UniformOutput', false ...
%                     )' ...
%                );
    
%     for firm=1:length(firms)
%         % Index for each customer of the firm.
%         idx = find(minInd==firm);
%         % Calculate the mean xy-coordinates / centroid of each firms market
%         centroid(firm,:) = mean(customers(idx,:),1);
%         % Calculate the distance form each customer to the respective market centroid.
%         distance_c(idx) = pdist2(centroid(firm,:), customers(idx,:), 'euclidean');
%     end

    %% Shape output variables
    l = sqrt(length(customers));
    customers = reshape(minInd,[l l]);
    utility = reshape(minVal,[l l]);
%     distance_centroid = reshape(distance_c,[l l]);
    
    % Output variables
    M = customers;
    U = utility;
%     C = centroid;
%     D = distance_centroid;
end