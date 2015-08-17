% Network -- simple model
% version 0.01
% Jonas K. Sekaman

clearvars;

% Preferences
pref.N = 3; % Number of firms
pref.iterations = 120;
pref.nv = 64; % Number of verticies/links between customers
pref.p = 0.2; % Rewiring probability in Erdos-Renyi network

%set(0,'DefaultFigureWindowStyle','docked')

% Setup
firm = 1:pref.N; % Gives each firm an ID/name
%rng('default'); % Seed such that the randomly generated results are repeatable

% Create Line network:
%network.Adj = zeros(pref.nv); network.Adj(logical(eye(size(network.Adj)))) = 1; network.Adj = tril(circshift(network.Adj, [1 0])); network.Adj = sparse(network.Adj); network.nv = pref.nv-1;
% Create Circle network:
%network.Adj = zeros(pref.nv); network.Adj(logical(eye(size(network.Adj)))) = 1; network.Adj = circshift(network.Adj, [1 0]); network.Adj = sparse(network.Adj); network.nv = pref.nv-1;
% Crate Erdos-Renyi network:
network = erdosRenyi(pref.nv, pref.p, 1); 

network.Adj = tril(network.Adj + network.Adj')~=0; % Convert to undirected network, where all non-zero elements equal 1.
network.Adj(logical(eye(size(network.Adj)))) = 0; % Remove any self-looping nodes.
symmetric = network.Adj + network.Adj';
ShortDistance = graphallshortestpaths(network.Adj, 'Directed', 'false');
ndegree = full(sum(symmetric,1)); % node degree
a = 0.3; % Fake alpha/transperency
color = linspecer(pref.N); % Generates a color for each of the firms that is evenly space out in / distinguishable.
colorAlpha = a*color + (1-a)*ones(pref.N,3); % Setting false alpha/transperency as weigheted average of color and white (since background is white)
% Initial position/node of firm
% TODO: Can the firms initially overlap, ie. have the same position? (change placement from true to false)
FirmE0 = randsample(1:network.nv, pref.N, true);

% Creating a 2D matrix of firms. 1st dimension is iteration, 2nd is firm, and value is the node position.
FirmE(1,:) = FirmE0;

% Matrix with each customer and their closest firm.
% Creating a 2D matrix of shares. 1st dimension is iteration, 2nd is node number, and value is closest firm.
share = marketsharen(FirmE(1,:), network.Adj);
shares(1,:) = share;
shares_mean(:,:,1) = histc(share(:), firm)/numel(share);
shares_degree(1,:) = marketdegree(share, symmetric);
%shares_centroid_distance

% Init image
h = view(biograph(network.Adj, [],'ShowArrows','off', 'LayoutType', 'equilibrium'));
set(h.Nodes, 'LineColor', [1 1 1]); % Change node line color to white
set(h.Edges,'LineColor',[0.9 0.9 0.9])
for j = 1:length(h.Nodes)
    h.Nodes(j).Label =...
			sprintf('%s',num2str(j)); % Set label equal to the node number

    if isempty(find(FirmE(1,:)==j))
        % if no firm on node
        h.Nodes(j).Color = colorAlpha(shares(1,j),:);
    else
        % if firm located on node
        h.Nodes(j).Color = color(shares(1,j),:);
        h.Nodes(j).FontSize = 12;
    end
end
h.ShowTextInNodes = 'label';
dolayout(h);


% Evolution
for i = 2:pref.iterations
    for n = 1:pref.N
        % Generate moves for firm n.
        moves = find(symmetric(FirmE(i-1,n),:));

        % Calculate share for each of the valid moves.
        clearvars shareSum;
        for m = 1:length(moves)
            FirmETrial = FirmE(i-1,:);
            % Replacing firm n's position with the position of move m.
            FirmETrial(n) = moves(m);
            % Calculating share with move, when other firms stay fixed.
            shareTrial = marketsharen(FirmETrial,network.Adj);
            shareSum(m) = mean(shareTrial == n);
        end
        % Find the move that gives the highest share
        [maxShare, maxInd] = max(shareSum);

        % Compare this with the current share
        if maxShare < shares_mean(n,:,i-1)
            % Stay
            FirmE(i,n) = FirmE(i-1,n);
        elseif maxShare > shares_mean(n,:,i-1)
            % Move
            FirmE(i,n) = moves(maxInd);
        else
            % Tie. Randomly draw
            FirmE(i,n) = cell2mat(randsample({FirmE(i-1,n), moves(maxInd)},1));
        end

    end
    share = marketsharen(FirmE(i,:),network.Adj);
    shares(i,:) = share;
    shares_mean(:,:,i) = histc(share(:), firm)/numel(share);
    shares_degree(i,:) = marketdegree(share, symmetric);
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
    a = ShortDistance(FirmE(i,:),FirmE(i,:));
    distance(i,:) = a(find(~tril(ones(size(a))))); % distance(i,:) = [12 13 14 23 24 34]
end
% Draw distaneces to other firms
figure(5);
maxd = max(distance(:));
clf reset; % Reset figure.
idx = 1:nchoosek(pref.N,2); % Creating an index array with all the combinations
for n = 1:pref.N-1 % Loop through each firm except the last one (since last firm  is only one that is present in all plots)
    subplot(pref.N-1,1,n); % Create stacked subplots
    set(gca, 'ColorOrder', color((n+1):pref.N,:), 'NextPlot', 'replacechildren'); % http://www.mathworks.com/matlabcentral/answers/19815#answer_26172
    %plot( distance(:,idx(1:pref.N-n)) ); % Using the index to plot repectively [12 13 14], [23 24], [34]
    stairs( distance(:,idx(1:pref.N-n)) ); % Using the index to plot repectively [12 13 14], [23 24], [34]
    set(gca,'ytick',0:maxd);
    set(gca,'ylim',[0,(maxd+1)]);
    legend(num2str([(n+1):pref.N]')); % Add legend with other firms numbers
    title(sprintf('Distance from firm %d to other firms',n)); % Add title
    idx(1:pref.N-n) = []; % Removing from the index the firms already plotted
end


% Draw firm position
%figure(6);
%clf reset; % Reset figure.
%hold on;
%title('Position of firms'); % Add title
%for n = 1:pref.N
%    stairs(FirmE(:,n), 'Color', color(n,:));
%    set(gca,'YTick',1:pref.nv);
%    set(gca,'YLim', [0,pref.nv+1]);
%end

figure(7);
[ndegree_s, ndegree_order] = sort(ndegree, 'descend');
[ndegree_sr, ndegree_order_r] = sort(ndegree_order);
clf reset; % Reset figure.
hold on;
title('Position of firms'); % Add title
for n = 1:pref.N
    stairs3(ndegree_order(FirmE(:,n)), 'Alpha', 0.05, 'Color', color(n,:));
    %stairs(ndegree_order(FirmE(:,n)), 'Color', color(n,:));
    set(gca,'YTick',1:pref.nv);
    set(gca,'YLim', [0,pref.nv+1]);
    set(gca,'YTickLabel',ndegree_order_r);
    set(gca,'YDir','reverse')
end
[C,ia,ic] = unique(ndegree_s, 'first');
for k = 1:length(ia)
    plot(repmat(ia(k), 1, pref.iterations), 'Color', [0.75 0.75 0.75], 'LineStyle', ':');
    text(pref.iterations, ia(k)+pref.nv*0.03, sprintf('Node degree %d',ndegree_s(ia(k))), 'Color', [0.75 0.75 0.75], 'HorizontalAlignment', 'right');
end
hold off;


% Degree of market
figure(8);
clf reset; % Reset figure.
hold on;
title('Degree of the market of the firm'); % Add title
for n = 1:pref.N
    stairs(shares_degree(:,n), 'Color', color(n,:));
end


% Draw evolution
h1 = view(biograph(network.Adj, [],'ShowArrows','off', 'LayoutType', 'equilibrium'));
set(h1.Edges,'LineColor',[0.9 0.9 0.9])
for j = 1:length(h1.Nodes)
    h1.Nodes(j).Label =...
			sprintf('%s',num2str(j)); % Set label equal to the node number
end
h1.ShowTextInNodes = 'label';
for i = 1:pref.iterations
    for j = 1:length(h1.Nodes)
        if isempty(find(FirmE(i,:)==j))
            % if no firm on node
            h1.Nodes(j).Color = colorAlpha(shares(i,j),:);
            h1.Nodes(j).LineColor = colorAlpha(shares(i,j),:);
            h1.Nodes(j).FontSize = 8;
        else
            % if firm located on node
            h1.Nodes(j).Color = color(shares(i,j),:);
            h1.Nodes(j).LineColor = color(shares(i,j),:);
            h1.Nodes(j).FontSize = 12;
        end
    end
    dolayout(h1, 'Pathsonly', true);
    pause(.1);
end