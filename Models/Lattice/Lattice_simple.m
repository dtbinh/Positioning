% Lattice -- simple model
% version 0.01
% Jonas K. Sekamane

% Inspired in part by: 
%   Ottino, B., Stonedahl, F. and Wilensky, U. (2009). 
%   NetLogo Hotelling?s Law model. 
%   http://ccl.northwestern.edu/netlogo/models/Hotelling'sLaw. 
%   Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

clearvars;

% Preferences
pref.length = 12; % Number of consumers on line, otherwise length of lattice on square
pref.height = 12; % Equal 1 if line, otherwise height of lattice on square
pref.N = 4; % Number of firms
pref.iterations = 400;

% Setup
firm = 1:pref.N; % Gives each firm an ID/name
rng('default'); % Seed such that the randomly generated results are repeatable
lattice = zeros(pref.height,pref.length);
% TODO: Can the firms initially overlap, ie. have the same position?
x0 = randi([1 pref.length], 1, pref.N); % Initial x-position of firm
y0 = randi([1 pref.height], 1, pref.N); % Initial y-position of firm
color = linspecer(pref.N); % Generates a color for each of the firms that is evenly space out in / distinguishable.

% Creating a 3D matrix of firms. 1st dimension is firm, 2nd is (x,y) coordinates, and 3rd is iteration.
xy = zeros(pref.N, 2, 2); 
xy(:,:,1) = [x0' y0'];

% Matrix with each customer and their closest firm.
% Creating a 3D matrix of shares. 1st dimension is x, 2nd is y, and 3rd is iteration.
shares = zeros(pref.height, pref.length, 2); 
share = marketshare(xy(:,:,1),lattice);
shares(:,:,1) = share;
shares_mean(:,:,1) = histc(share(:), firm)/numel(share);
share_props = regionpropsext(share, 'Centroid', 'Extrema', 'Perimeter');
shares_perimeter(:,:,1) = cat(1,share_props.Perimeter);
shares_extrema_perimeter(:,:,1) = cat(1,share_props.ExtremaPerimeter);
shares_centroid_distance(:,:,1) = diag(pdist2(xy(:,:,1), cat(1,share_props.Centroid), 'Euclidean'));


% Init image
figure(1);
clf reset; % Reset figure.
image(shares(:,:,1),'AlphaData', 0.2); % Create image of the marketshare with dimmed opacity.
title('Market initially'); % Add title
colormap(color); % Set the respective firm colors on each market share.
axis image; % Scale to fit the dimensions of customers/lattice.
set(gca,'ytick', 1:pref.height); % Adjust axis to integer
set(gca,'xtick', 1:pref.length); % Alternative change 1:pref.length to [] to remove axis labels.
set(gca,'ydir', 'normal'); % By default image() reverese the y-axis. This reinstates the normal ordering.
hold on; % Place following scatter/plots on top of image.
scatter( xy(:,1,1,:)' , xy(:,2,1,:)' , [], color, 'filled'); % Plot the firms with respective colors.
hold off; % Don't place anymore plots on top of figure.


% Evolution
for i = 2:pref.iterations
   for n = 1:pref.N
        % Generate moves for firm n.
        N = xy(n,:,i-1,:) + [0 1];
        E = xy(n,:,i-1,:) + [1 0];
        S = xy(n,:,i-1,:) + [0 -1];
        W = xy(n,:,i-1,:) + [-1 0];

        moves = [N' E' S' W']; % Combining the moves. First row is new x-values, second is y-values

        % Check that moves are in range.
        % Y-axis in  range
        idx = find(moves(1,:) >= 1 & moves(1,:) <= pref.length, size(moves,2) );
        moves = moves(:,idx);
        % X-axis in range
        idx = find(moves(2,:) >= 1 & moves(2,:) <= pref.height, size(moves,2) );
        moves = moves(:,idx);

        % Calculate share for each of the remaing valid moves.
        clearvars shareSum;
        for m = 1:size(moves,2)
            xyTrial = xy(:,:,i-1);
            % Replacing firm n's coordinates with the coordinates of move m.
            xyTrial(n,:) = moves(:,m)';
            % Calculating share with move, when other firms stay fixed.
            shareTrial = marketshare(xyTrial,lattice);
            shareSum(m) = mean(shareTrial(:) == n);
        end
        % Find the move that gives the highest share
        [maxShare, maxInd] = max(shareSum);
        
        % Compare this with the current share
        if maxShare < shares_mean(n,:,i-1)
            % Stay
            xy(n,:,i,:) = xy(n,:,i-1,:);
        elseif maxShare > shares_mean(n,:,i-1)
            % Move
            xy(n,:,i,:) = moves(:,maxInd);
        else
            % Tie. Randomly draw
            xy(n,:,i,:) = cell2mat(randsample({xy(n,:,i-1,:)', moves(:,maxInd)},1));
        end
        
   end
   
   share = marketshare(xy(:,:,i),lattice);
   shares(:,:,i) = share;
   shares_mean(:,:,i) = histc(share(:), firm)/numel(share);
   share_props = regionpropsext(share, 'Centroid', 'Extrema', 'Perimeter');
   shares_perimeter(:,:,i) = cat(1,share_props.Perimeter);
   shares_extrema_perimeter(:,:,i) = cat(1,share_props.ExtremaPerimeter);
   shares_centroid_distance(:,:,i) = diag(pdist2(xy(:,:,i), cat(1,share_props.Centroid), 'Euclidean'));
   
end


% Draw market share
figure(3);
clf reset; % Reset figure.
hold on;
title('Share of market'); % Add title
for n = 1:pref.N
    plot(squeeze(shares_mean(n,1,:)), 'Color', color(n,:));
end
plot(repmat(1/pref.N, 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');

figure(4);
clf reset; % Reset figure.
hold on;
title('Share of market moving average (lagged 10% of number of iterations)'); % Add title
for n = 1:pref.N
    simple = tsmovavg(squeeze(shares_mean(n,1,:)), 's', pref.iterations*0.1, 1);
    plot(simple, 'Color', color(n,:), 'LineStyle', '--');
end
plot(repmat(1/pref.N, 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');


% Calculate the distance to each of the other firms
for i = 1:pref.iterations
    distance(i,:) = pdist(xy(:,:,i), 'euclidean'); % distance(i,:) = [12 13 14 23 24 34]
end
% Draw distaneces to other firms
figure(5);
clf reset; % Reset figure.
idx = 1:nchoosek(pref.N,2); % Creating an index array with all the combinations
for n = 1:pref.N-1 % Loop through each firm except the last one (since last firm  is only one that is present in all plots)
    subplot(pref.N-1,1,n); % Create stacked subplots
    set(gca, 'ColorOrder', color((n+1):pref.N,:), 'NextPlot', 'replacechildren'); % http://www.mathworks.com/matlabcentral/answers/19815#answer_26172
    plot( distance(:,idx(1:pref.N-n)) ); % Using the index to plot repectively [12 13 14], [23 24], [34]
    legend(num2str([(n+1):pref.N]')); % Add legend with other firms numbers
    title(sprintf('Distance from firm %d to other firms',n)); % Add title
    idx(1:pref.N-n) = []; % Removing from the index the firms already plotted
end

% Draw distance to centroids
figure(6);
clf reset; % Reset figure.
hold on;
title('Firm distance to the firm`s market centroid'); % Add title
for n = 1:pref.N
    plot(squeeze(shares_centroid_distance(n,1,:)), 'Color', color(n,:));
end

% Draw perimeter
%figure(7);
%clf reset; % Reset figure.
%hold on;
%title('Market Perimeter'); % Add title
%for n = 1:pref.N
%    plot(squeeze(shares_perimeter(n,1,:)), 'Color', color(n,:));
%end

figure(8);
clf reset; % Reset figure.
hold on;
title('Perimeter of market extrema points'); % Add title
for n = 1:pref.N
    plot(squeeze(shares_extrema_perimeter(n,1,:)), 'Color', color(n,:));
end

% Equal split market
%figure(9);
%clf reset; % Reset figure.
%hold on;
%title('Equal split market'); % Add title
%for n = 1:pref.N
%    plot(abs(squeeze(shares_mean(3,1,:))-1/pref.N).*squeeze(shares_centroid_distance(n,1,:)), 'Color', color(n,:));
%end



% Draw evolution
for i = 1:pref.iterations
    figure(2);
    clf reset; % Reset figure.
    image(shares(:,:,i),'AlphaData', 0.2); % Create image of the marketshare with dimmed opacity.
    title(sprintf('Market (iteration %d)',i)); % Add title
    colormap(color); % Set the respective firm colors on each market share.
    axis image; % Scale to fit the dimensions of customers/lattice.
    set(gca,'ytick', []); % Adjust axis to integer
    set(gca,'xtick', []); % Alternative change 1:pref.length to [] to remove axis labels.
    set(gca,'ydir', 'normal'); % By default image() reverese the y-axis. This reinstates the normal ordering.
    hold on; % Place following scatter/plots on top of image.
    scatter( xy(:,1,i,:)' , xy(:,2,i,:)' , [], color, 'filled'); % Plot the firms with respective colors.
    for n = 1:pref.N
       n_xy = xy(n,:,1:i,:);
       n_xy = reshape(n_xy, [size(n_xy,2) size(n_xy,3) 1]);
       line(n_xy(1,:)', n_xy(2,:)', 'Color', color(n,:), 'LineStyle', ':');
    end
    hold off; % Don't place anymore plots on top of figure.
    pause(.1);
end
