function F = marketdegree(share, symmetric)
%marketdegree   calculates the degree (links going in/out) of the market for each firm
%   MARKETDEGREE(SHARE, SYMMETRIC) 
%   where SHARE is a list of all nodes with the value indicating the firm.
%   and SYMMETRIC is a symmetric undirected adjecncy matrix
%   Jonas K. Sekamane. 
%   Version 0.01

    Max = max(share(:));
    % Array with each customer and their closest firm.
    Degree = [];
    for n = 1:Max
        A = symmetric(share==n,:);
        B = A(:,share~=n);
        Degree = [Degree sum(B(:))]; % Alternatively: sum(sum(B,1)~=0)
    end

    F = full(Degree);
end