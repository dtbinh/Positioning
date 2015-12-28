%% PLOTS
% Jonas K. Sekamane
%
% Various plots used to vertify partical results and debug code

% Draw evolution
for i = 1:pref.iterations
    figure(2);
    clf reset; % Reset figure.
    image(x, y, market(:,:,i),'AlphaData', 0.2); % Create image of the marketshare with dimmed opacity.
    title(sprintf('Market (iteration %d)',i)); % Add title
    colormap(color); % Set the respective firm colors on each market share.
    axis image; % Scale to fit the dimensions of customers/lattice.
    set(gca,'ydir', 'normal'); % By default image() reverese the y-axis. This reinstates the normal ordering.
    hold on; % Place following scatter/plots on top of image.
    scatter( xy(:,1,i,:)' , xy(:,2,i,:)' , [], color, 'filled'); % Plot the firms with respective colors.
%    scatter(centroid(:,1,i)',centroid(:,2,i)', [], color, 'x'); % Plot current market centroid point
%     for n = 1:pref.N % Trace the firm movement over time
%        n_xy = xy(n,:,1:i,:);
%        n_xy = reshape(n_xy, [size(n_xy,2) size(n_xy,3) 1]);
%        line(n_xy(1,:)', n_xy(2,:)', 'Color', color(n,:), 'LineStyle', ':');
%     end
    hold off; % Don't place anymore plots on top of figure.
    pause(.1);
end


% 3D plot of population density function
figure(10);
surf(x,y,F, 'EdgeColor', 'none');%,'EdgeColor','none');
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-5 5 -5 5 0 .5])
title('Population density');
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');


% 2D plot of contours of population and firm initial position
figure(11);
contour(x,y,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
hold on;
scatter( xy_0(:,1)' , xy_0(:,2)', [], color, 'filled'); % Plot the firms with respective colors.
xlim([-5 5]); ylim([-5 5]);
title('Population contours and initial firm position');
xlabel('x1'); ylabel('x2');
hold off;


% Voronoi diagram of initial market
figure(12);
voronoi(xy_0(:,1)' , xy_0(:,2)');
xlim([-pref.boundary/2 pref.boundary/2]); ylim([-pref.boundary/2 pref.boundary/2]);
hold on;
scatter(centroid(:,1,1)',centroid(:,2,1)', [], color, 'x');
title('Voronoi diagram of initial market (with centroids)');
hold off;


% Plot of customer utility
figure(14);
imagesc(utility_i)
set(gca,'ydir', 'normal');
cm=flipud(jet);
colormap(cm);
title('Customer (dis)utility at final interation');


% Final market
figure(13);
h = imagesc(x,y, market_i);
set(h,'AlphaData',0.5);
set(gca,'ydir', 'normal');
colormap(color); % Set the respective firm colors on each market share.
hold on;
scatter( centroid(:,1,pref.iterations)' , centroid(:,2,pref.iterations)', [], color, 'x');
scatter( xy(:,1,pref.iterations,:)' , xy(:,2,pref.iterations,:)' , [], color, 'filled');
hold off;
title('Final market');


% Compass of the direction of the firm's movement
% 1. Trace past directions:
% figure(15);
% coloralpha = [color repmat(1/pref.iterations, 1, pref.N)'];
% hold on;
% for i = 1:pref.iterations
%     [compass_x, compass_y] = pol2cart(heading(i,:), ones(1,5));
%     h = compass(compass_x, compass_y);
%     set(h, {'Color'}, num2cell(coloralpha,2), 'LineWidth',2)
%     pause(.1);
% end
% hold off;
% 2. Show the current direction:
for i = 1:pref.iterations
    figure(15);
    [compass_x, compass_y] = pol2cart(heading(i,:), ones(1,pref.N));
    h = compass(compass_x, compass_y);
    set(h, {'Color'}, num2cell(color,2), 'LineWidth',2);
    title(sprintf('Heading (iteration %d)',i)); % Add title
    pause(.1);
end


% Draw market share
figure(3);
clf reset; % Reset figure.
hold on;
title('Share of market'); % Add title
for n = 1:pref.N
    plot(shares(:,n,:), 'Color', color(n,:));
end
plot(repmat(1/pref.N, 1, pref.iterations), 'Color', 'k', 'LineStyle', ':'); % Equal shares
hold off;

figure(4);
clf reset; % Reset figure.
hold on;
title('Share of market moving average (lagged approx. 10% of number of iterations)'); % Add title
for n = 1:pref.N
    simple = tsmovavg(shares(:,n,:), 's', round(pref.iterations*0.1), 1);
    plot(simple, 'Color', color(n,:));
end
plot(repmat(1/pref.N, 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');
hold off;


% Draw distance to centroids
figure(6);
clf reset; % Reset figure.
hold on;
title('Firm distance to the firm`s market centroid'); % Add title
for n = 1:pref.N
    plot(centroid_distance(:,n,:), 'Color', color(n,:));
end
hold off;


% Mean eccentricity
figure(16);
clf reset; % Reset figure.
plot(mean_eccentricity);
title('Mean eccentricity'); % Add title


% Effective Number of Parties/Players/Firms (ENP)
figure(17);
clf reset; % Reset figure.
plot(ENP);
title('ENP'); % Add title
