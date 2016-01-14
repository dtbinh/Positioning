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
%    scatter(centroid(:,1,i)',centroid(:,2,i)', [], color, 'x'); % Plot current market centroid point
    for n = 1:pref.N % Trace the firm movement over time
       n_xy = xy(n,:,1:i,:);
       n_xy = reshape(n_xy, [size(n_xy,2) size(n_xy,3) 1]);
       line(n_xy(1,:)', n_xy(2,:)', 'Color', color(n,:), 'LineStyle', ':');
    end
    scatter( xy(:,1,i,:)' , xy(:,2,i,:)' , [], color, 'filled'); % Plot the firms with respective colors.
    hold off; % Don't place anymore plots on top of figure.
    pause(.03);
end


% 3D plot of population density function
figure(10);
surf(x,y,F);%,'EdgeColor','none');
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-5 5 -5 5 0 .2])
%title('Population density');
title({'Population density', ['\mu_r = -\mu_l = ', num2str(pref.mu) ';   n_l/n_r: ' num2str(pref.n_ratio)]}); % Add title
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
set(gca,'ZTick',[0:0.05:0.2]);


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

% Plot of customer utility
figure(141);
imagesc(utility_i.*F)
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
    pause(.02);
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
hold on;
m2nd_mean_eccentricity = mean(mean_eccentricity(length(mean_eccentricity)/2+1:end));
sd2nd_mean_eccentricity = std(mean_eccentricity(length(mean_eccentricity)/2+1:end));
plot(repmat( m2nd_mean_eccentricity , 1, pref.iterations), 'k');
plot(repmat( m2nd_mean_eccentricity + sd2nd_mean_eccentricity , 1, pref.iterations), 'k:');
plot(repmat( m2nd_mean_eccentricity - sd2nd_mean_eccentricity , 1, pref.iterations), 'k:');
hold off;


% Effective Number of Parties/Players/Firms (ENP)
figure(17);
clf reset; % Reset figure.
plot(ENP);
title('ENP'); % Add title
hold on;
m2nd_ENP = mean(ENP(length(mean_eccentricity)/2:end));
sd2nd_ENP = std(ENP(length(mean_eccentricity)/2:end));
plot(repmat( m2nd_ENP , 1, pref.iterations), 'k');
plot(repmat( m2nd_ENP + sd2nd_ENP , 1, pref.iterations), 'k:');
plot(repmat( m2nd_ENP - sd2nd_ENP , 1, pref.iterations), 'k:');
hold off;


% Mean eccentricity
figure(18);
clf reset; % Reset figure.
plot(mean_eccentricity_1(3,:)');
title('Mean eccentricity'); % Add title
hold on;
m2nd_mean_eccentricity = mean(mean_eccentricity(length(mean_eccentricity)/2+1:end));
sd2nd_mean_eccentricity = std(mean_eccentricity(length(mean_eccentricity)/2+1:end));
plot(repmat( m2nd_mean_eccentricity , 1, pref.iterations), 'k');
plot(repmat( m2nd_mean_eccentricity + sd2nd_mean_eccentricity , 1, pref.iterations), 'k:');
plot(repmat( m2nd_mean_eccentricity - sd2nd_mean_eccentricity , 1, pref.iterations), 'k:');
hold off;

% Effective Number of Parties/Players/Firms (ENP)
figure(180);
clf reset; % Reset figure.
plot(ENP);
title('ENP'); % Add title
hold on;
m2nd_ENP = mean(ENP(length(ENP)/2:end));
sd2nd_ENP = std(ENP(length(ENP)/2:end));
plot(repmat( m2nd_ENP , 1, pref.iterations), 'k');
plot(repmat( m2nd_ENP + sd2nd_ENP , 1, pref.iterations), 'k:');
plot(repmat( m2nd_ENP - sd2nd_ENP , 1, pref.iterations), 'k:');
hold off;

% Effective Number of Parties/Players/Firms (ENP)
figure(181);
clf reset; % Reset figure.
plot(ENP_1(13,:)');
title('ENP'); % Add title
hold on;
m2nd_ENP = mean(ENP(length(ENP)/2:end));
sd2nd_ENP = std(ENP(length(ENP)/2:end));
plot(repmat( m2nd_ENP , 1, pref.iterations), 'k');
plot(repmat( m2nd_ENP + sd2nd_ENP , 1, pref.iterations), 'k:');
plot(repmat( m2nd_ENP - sd2nd_ENP , 1, pref.iterations), 'k:');
hold off;

% Mean eccentricity kdensity
figure(19);
for rep = 1:pref.repetitions
    [f,xi] = ksdensity(mean_eccentricity_1(rep,:));
    plot(xi,f);
    hold on;
end
title('Mean eccentricity'); % Add title
hold off;

figure(20);
for rep = 1:pref.repetitions
    [f,xi] = ksdensity(ENP_1(rep,:));
    plot(xi,f);
    hold on;
end
title('ENP'); % Add title
hold off;

figure(21);
[f,xi] = ksdensity(ENP);
plot(xi,f);
title('ENP'); % Add title

figure(22);
[f,xi] = ksdensity(mean_eccentricity);
plot(xi,f);
title('Mean eccentricity'); % Add title

figure(210);
[f,xi] = ksdensity(ENP_1(3,:));
plot(xi,f);
title('ENP'); % Add title


run = 6;
rep = 3;

figure(23);
plot(data_ENP_within(rep,:,run)');
title('ENP within +/- 1 std. dev.'); % Add title

figure(24);
plot(data_ENP(rep,:,run)');
hold on;
plot(repmat( est_ENP(rep,:,run) , 1, pref.iterations), 'k');
plot(repmat( est_ENP(rep,:,run) + est_std_ENP(rep,:,run) , 1, pref.iterations), 'k:');
plot(repmat( est_ENP(rep,:,run) - est_std_ENP(rep,:,run) , 1, pref.iterations), 'k:');
title('ENP'); % Add title
hold off;

figure(241);
plot(data_mean_eccentricity(rep,:,run)');
hold on;
plot(repmat( est_mean_eccentricity(rep,:,run) , 1, pref.iterations), 'k');
plot(repmat( est_mean_eccentricity(rep,:,run) + est_std_mean_eccentricity(rep,:,run) , 1, pref.iterations), 'k:');
plot(repmat( est_mean_eccentricity(rep,:,run) - est_std_mean_eccentricity(rep,:,run) , 1, pref.iterations), 'k:');
title('Mean eccentricity'); % Add title
hold off;




%% Plot mixed
% High contrast colours:
%ix = randperm(64); % Matlab colors contain 64 values
%jettemp = jet;
%jetshuffle = jettemp(ix,:);



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
%    scatter(centroid(:,1,i)',centroid(:,2,i)', [], color, 'x'); % Plot current market centroid point
%     for n = 1:pref.N % Trace the firm movement over time
%        n_xy = xy(n,:,1:i,:);
%        n_xy = reshape(n_xy, [size(n_xy,2) size(n_xy,3) 1]);
%        line(n_xy(1,:)', n_xy(2,:)', 'Color', color(n,:), 'LineStyle', ':');
%     end
    scatter( xy(:,1,i,:)' , xy(:,2,i,:)', [], colorrulefirm, 'filled'); % Plot the firms with respective colors.
    hold off; % Don't place anymore plots on top of figure.
    pause(.03);
end


% Effective Number of Parties/Players/Firms (ENP)
figure(17);
clf reset; % Reset figure.
plot(ENP);
title('ENP'); % Add title
hold on;
m2nd_ENP = mean(ENP(length(mean_eccentricity)/2:end));
sd2nd_ENP = std(ENP(length(mean_eccentricity)/2:end));
plot(repmat( m2nd_ENP , 1, pref.iterations/pref.psi), 'k');
plot(repmat( m2nd_ENP + sd2nd_ENP , 1, pref.iterations/pref.psi), 'k:');
plot(repmat( m2nd_ENP - sd2nd_ENP , 1, pref.iterations/pref.psi), 'k:');
hold off;


% Mean eccentricity
figure(16);
clf reset; % Reset figure.
plot(mean_eccentricity, 'k');
title('Mean eccentricity'); % Add title
hold on;
for r=1:length(pref.ruleset)
    plot(mean_eccentricity_rule(:,r), 'Color', colorrule(r,:));
end
m2nd_mean_eccentricity = mean(mean_eccentricity(length(mean_eccentricity)/2+1:end));
sd2nd_mean_eccentricity = std(mean_eccentricity(length(mean_eccentricity)/2+1:end));
plot(repmat( m2nd_mean_eccentricity , 1, pref.iterations/pref.psi), 'k');
plot(repmat( m2nd_mean_eccentricity + sd2nd_mean_eccentricity , 1, pref.iterations/pref.psi), 'k:');
plot(repmat( m2nd_mean_eccentricity - sd2nd_mean_eccentricity , 1, pref.iterations/pref.psi), 'k:');
hold off;

% Mean share
figure(161);
clf reset; % Reset figure.
plot(mean_share, 'k');
title('Mean market share'); % Add title
hold on;
for r=1:length(pref.ruleset)
    plot(mean_share_rule(:,r), 'Color', colorrule(r,:));
end
m2nd_mean_share = mean(mean_share(length(mean_share)/2+1:end));
sd2nd_mean_share = std(mean_share(length(mean_share)/2+1:end));
plot(repmat( m2nd_mean_share , 1, pref.iterations/pref.psi), 'k');
plot(repmat( m2nd_mean_share + sd2nd_mean_share , 1, pref.iterations/pref.psi), 'k:');
plot(repmat( m2nd_mean_share - sd2nd_mean_share , 1, pref.iterations/pref.psi), 'k:');
hold off;


% Number of surviving firms
figure(171);
clf reset; % Reset figure.
plot(N,'k');
title('Number of surviving firms'); % Add title
hold on;
for r=1:length(pref.ruleset)
    plot(N_rule(:,r), 'Color', colorrule(r,:));
end
m2nd_N = mean(N(length(N)/2:end));
sd2nd_N = std(N(length(N)/2:end));
plot(repmat( m2nd_N , 1, pref.iterations/pref.psi), 'k');
plot(repmat( m2nd_N + sd2nd_N , 1, pref.iterations/pref.psi), 'k:');
plot(repmat( m2nd_N - sd2nd_N , 1, pref.iterations/pref.psi), 'k:');
hold off;


% Mean age at death
figure(162);
clf reset; % Reset figure.
plot(mean_age_death,'k');
title('Mean firm age at death'); % Add title
hold on;
for r=1:length(pref.ruleset)
    plot(mean_age_death_rule(:,r), 'Color', colorrule(r,:));
end
%m_age = mean(mean_age_death, 'omitnan');
%plot(repmat( m_age , 1, pref.iterations/pref.psi), 'k');
%sd2nd_mean_share = std(mean_share(length(mean_share)/2+1:end));
%plot(repmat( m2nd_mean_share + sd2nd_mean_share , 1, pref.iterations/pref.psi), 'k:');
%plot(repmat( m2nd_mean_share - sd2nd_mean_share , 1, pref.iterations/pref.psi), 'k:');
hold off;

