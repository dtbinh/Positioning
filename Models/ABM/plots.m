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
    pause(.01);
end


% 3D plot of population density function
figure(10);
surf(x,y,F);%,'EdgeColor','none');
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-3 3 -3 3 0 .2])
%title('Population density');
title({'Population density', ['\mu_r = -\mu_l = ', num2str(pref.mu) ';   n_l/n_r: ' num2str(pref.n_ratio)]}); % Add title
xlabel('x'); ylabel('y'); zlabel('Probability density');
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

% 2D plot of contours of population and firm initial position
figure(11);
contour(x,y,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
hold on;
scatter( xy(:,1,pref.iterations,:)' , xy(:,2,pref.iterations,:)', [], color, 'filled'); % Plot the firms with respective colors.
xlim([-5 5]); ylim([-5 5]);
title('Population contours and firm position');
xlabel('x'); ylabel('y');
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
title('Market areas, market centroids and firm position');
xlabel('x'); ylabel('y');


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
%plot(repmat( m2nd_ENP , 1, pref.iterations), 'k');
%plot(repmat( m2nd_ENP + sd2nd_ENP , 1, pref.iterations), 'k:');
%plot(repmat( m2nd_ENP - sd2nd_ENP , 1, pref.iterations), 'k:');
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
for i = 1:pref.iterations*pref.psi
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
plot(repmat( m2nd_ENP , 1, pref.iterations), 'k');
plot(repmat( m2nd_ENP + sd2nd_ENP , 1, pref.iterations), 'k:');
plot(repmat( m2nd_ENP - sd2nd_ENP , 1, pref.iterations), 'k:');
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
plot(repmat( m2nd_mean_eccentricity , 1, pref.iterations), 'k');
plot(repmat( m2nd_mean_eccentricity + sd2nd_mean_eccentricity , 1, pref.iterations), 'k:');
plot(repmat( m2nd_mean_eccentricity - sd2nd_mean_eccentricity , 1, pref.iterations), 'k:');
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
plot(repmat( m2nd_mean_share , 1, pref.iterations), 'k');
plot(repmat( m2nd_mean_share + sd2nd_mean_share , 1, pref.iterations), 'k:');
plot(repmat( m2nd_mean_share - sd2nd_mean_share , 1, pref.iterations), 'k:');
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
plot(repmat( m2nd_N , 1, pref.iterations), 'k');
plot(repmat( m2nd_N + sd2nd_N , 1, pref.iterations), 'k:');
plot(repmat( m2nd_N - sd2nd_N , 1, pref.iterations), 'k:');
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



%% Plot Inductor

% Draw evolution
for i = 1:pref.iterations*pref.psi
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
    scatter( xy(:,1,i,:)' , xy(:,2,i,:)' , [], color, 'filled'); % Plot the firms with respective colors.
    hold off; % Don't place anymore plots on top of figure.
    xlim([-5 5]); ylim([-5 5]);
    pause(.01);
end

% Draw market share
figure(3);
clf reset; % Reset figure.
hold on;
title('Share of market'); % Add title
for n = 1:pref.N
    plot(shares(:,n,:), 'Color', color(n,:));
end
plot(repmat(1/pref.N, 1, pref.iterations*pref.psi), 'Color', 'k', 'LineStyle', ':'); % Equal shares
hold off;

figure(4);
clf reset; % Reset figure.
hold on;
title('Share of market moving average (lagged approx. 10% of number of iterations)'); % Add title
for n = 1:pref.N
    simple = tsmovavg(shares(:,n,:), 's', round(pref.iterations*pref.psi*0.1), 1);
    plot(simple, 'Color', color(n,:));
end
plot(repmat(1/pref.N, 1, pref.iterations*pref.psi), 'Color', 'k', 'LineStyle', ':');
hold off;

% Check accuracy
accuracy_final = squeeze(cf_i(:, 24, :));
accuracy_final(accuracy_final==0) = NaN; % Non-zero accuracy.

accuracy_i = squeeze(cf(:, 24, :, :));
accuracy_positive = accuracy_i;
accuracy_positive(accuracy_positive==0) = NaN; % Non-zero accuracy.

accuracy_positive_mean_iter = mean(accuracy_positive, 3, 'omitnan');  % Mean over all iterations.
accuracy_positive_mean_firm = squeeze(mean(accuracy_positive, 1, 'omitnan'));  % Mean over all firm's condition/forcast rules.

% Active
active_i = squeeze(cf(:, 25, :, :));
active_positive = active_i;
active_positive(accuracy_positive==0) = NaN; % Non-zero accuracy.
active_positive_mean_firm = squeeze(mean(active_positive, 1, 'omitnan'));  % Mean over all firm's condition/forcast rules.

% accuracy "per times active".
accuracy_active_positive = accuracy_positive./(1+active_positive);
accuracy_active_positive_mean_firm = squeeze(mean(accuracy_active_positive, 1, 'omitnan'));

% Accuracy excluding default rule
accuracyd_i = squeeze(cf(2:end, 24, :, :));
accuracyd_positive = accuracyd_i;
accuracyd_positive(accuracyd_positive==0) = NaN; % Non-zero accuracy.
accuracyd_positive_mean_firm = squeeze(mean(accuracyd_positive, 1, 'omitnan'));  % Mean over all firm's condition/forcast rules.


figure(901);
histogram(accuracy_final);
title('Accuracy at final iteration');

figure(902);
histogram(accuracy_positive_mean_iter);
title('Mean accuracy over all iteration');

figure(903);
clf reset; % Reset figure.
%plot(accuracy_positive_mean_firm');
hold on; for n=1:pref.N plot(accuracy_positive_mean_firm(n,1:500), 'Color', color(n,:)); end; hold off;
title('Mean accuracy of all the firm´s condition/forecast rules'); % Add title
%hold on;
%plot(mean(accuracy_positive_mean_firm,1), 'k');
%hold off;

figure(904);
clf reset; % Reset figure.
plot(mean(accuracy_positive_mean_firm,1), 'k');
title('Mean accuracy of all condition/forecast rules'); % Add title

figure(913);
clf reset; % Reset figure.
hold on;
plot(sqrt(mean_forecasterror(1:500,:)));
hold off;
title('Mean forecast error of all firm'); % Add title

figure(914);
clf reset; % Reset figure.
plot(accuracy_active_positive_mean_firm');
title('Mean accuracy/(1+active) of all the firm´s condition/forecast rules');

figure(915);
clf reset; % Reset figure.
plot(accuracyd_positive_mean_firm');
title('Mean accuracy of all the firm´s condition/forecast rules (excluding default)'); % Add title




% Which conditions are being used?

% % Total number of times that a state is active.
sum(~isnan(accuracy_final)) % Percentage of conditions that have been used.
figure(9001);
n=1;
stem(~isnan(accuracy_final(:,n)));
title('Total number of states occurrences'); % Add title

% % Total number of times that a state is active.
state_occurrences = squeeze(sum(J_states,1,'omitnan'));
state_occurrences_total = sum(state_occurrences,2,'omitnan');
figure(905);
stem(state_occurrences_total);
title('Total number of states occurrences'); % Add title

state_occurrences_type = [mean(state_occurrences(1:5,:)); mean(state_occurrences(6:9,:)); mean(state_occurrences(10:11,:)); state_occurrences(12,:); state_occurrences(13,:)];
figure(906);
plot(state_occurrences_type(:,1:500)');
title('State occurrences over time'); % Add title

size(state_occurrences_type)

% figure(9060);
% plot(state_occurrences(6:8,1:50)');
% title('MA-state occurrences over time'); % Add title
% figure(9061);
% plot(state_occurrences(9:11,1:50)');
% title('MA-state occurrences over time'); % Add title


% Total number of times that a state is active. second half
state2nd_occurrences = squeeze(sum(J_states(:,:,101:200),1,'omitnan'));
state2nd_occurrences_total = sum(state2nd_occurrences,2,'omitnan');
figure(907);
stem(state2nd_occurrences_total);
title('Total number of states occurrences'); % Add title

state2nd_occurrences_type = [mean(state2nd_occurrences(1:5,:)); mean(state2nd_occurrences(6:9,:)); mean(state2nd_occurrences(10:11,:)); state2nd_occurrences(12,:); state2nd_occurrences(13,:)];
figure(908);
plot(state2nd_occurrences_type');
title('State occurrences over time'); % Add title



size(J_states) % firm, bits, iteration
size(cf) % cf-rule, conditions/values, firms, iteration

bitsset = mean(~isnan(cf(:,1:13,:,:)), 1);
size(bitsset)

bitsset_funda = mean(bitsset(:,1:5,:,:), 2);
bitsset_trend = mean(bitsset(:,6:9,:,:), 2);
bitsset_ocsil = mean(bitsset(:,10:11,:,:), 2);
bitsset_allon = mean(bitsset(:,12,:,:), 2);
bitsset_aloff = mean(bitsset(:,13,:,:), 2);

% Faction of condition bits set to respectively fundamental, trending
% ocsillation and on/off.
figure(1001);
clf reset; hold on; for n=1:pref.N plot(squeeze(bitsset_funda(:,:,n,:)), 'Color', color(n,:)); end; hold off;
figure(1002);
clf reset; hold on; for n=1:pref.N plot(squeeze(bitsset_trend(:,:,n,:)), 'Color', color(n,:)); end; hold off;
figure(1003);
clf reset; hold on; for n=1:pref.N plot(squeeze(bitsset_ocsil(:,:,n,:)), 'Color', color(n,:)); end; hold off;
figure(1004);
plot(squeeze(bitsset_allon)');
figure(1005);
plot(squeeze(bitsset_aloff)');
% Decrease due to specificity cost.

cf_used_change = diff(cf(:,25,:,:), 1, 4);
cf_used_list = (cf_used_change > 0);

%size(cf_used_list)
figure(1011);
imagesc( squeeze(cf_used_list(:,:,1,1:500)) );
set(gca,'ydir', 'normal');

figure(1012);
imagesc( squeeze(cf_used_list(:,:,2,1:500)) );
set(gca,'ydir', 'normal');

figure(1013);
imagesc( squeeze(cf_used_list(:,:,3,1:500)) );
set(gca,'ydir', 'normal');

figure(1014);
imagesc( squeeze(cf_used_list(:,:,4,1:500)) );
set(gca,'ydir', 'normal');

figure(1015);
imagesc( squeeze(cf_used_list(:,:,5,1:500)) );
set(gca,'ydir', 'normal');



% If want to compare with maxcov-inductor-GA, then count the number of
% unique CF-rules used in each interval of psi interations. If used
% CF-rules converge to a limited set, then the number of unique rules
% should decrease.


% cf2 = cf;
% cf2(:,27,:,2:end) = cf_used_list;
% %size(cf_used_change)
% 
% cf_used = nan(1, size(cf2,2), size(cf2,3), size(cf2,4));
% for c=1:size(cf2,1)
%     for n=1:size(cf2,3)
%         for i=1:size(cf2,4)
%             if (cf2(c,27,n,i) == 1)
%               cf_used(1,:,n,i) = cf2(c,:,n,i);
%             end
%         end
%     end
% end

used_bitsset = ~isnan(cf_used(:,1:13,:,:));
used_bitsset_funda = mean(used_bitsset(:,1:5,:,:), 2);
used_bitsset_trend = mean(used_bitsset(:,6:9,:,:), 2);
used_bitsset_ocsil = mean(used_bitsset(:,10:11,:,:), 2);

figure(1011);
plot(squeeze(used_bitsset_funda)');
figure(1012);
plot(squeeze(used_bitsset_trend)');
figure(1013);
plot(squeeze(used_bitsset_ocsil)');


% % Find used cf
% I = find(cf_used_change > 0);
% dims = size(cf_used_change);
% smax = cell(size(dims));
% [smax{:}] = ind2sub(dims, I);
% %smax{1}
% %smax{4}
% dims2 = size(cf);
% dims2([1 3 4])
% idx_used = sub2ind(dims2([1 3 4]), smax{1}, smax{3}, smax{4});

cf_used_change = diff(cf(:,25,:,:), 1, 4);
cf_used_list = (cf_used_change > 0);
cf2 = cf;
cf2(:,27,:,2:end) = cf_used_list;
j = 1;
test_used = nan(1, 29);
for c=1:size(cf2,1)
    for n=1:size(cf2,3)
        for i=1:size(cf2,4)
            if (cf2(c,27,n,i) == 1)
              test_used(j,:) = [cf2(c,1:26,n,i) i n c];
              j = j+1;
            end
        end
    end
end

csvwrite(strcat('data/maxcov-inductor-ga_N5_mu12_n15_usedcf.csv'), test_used);



%cf_used2 = cf_used;

%delete_cf = zeros(size(cf_used,1), size(cf_used,3), size(cf_used,4));
%size(delete_cf)
%size(cf_used2)

%[~, idx] = sort(cf_used, 2);
%size(idx)
%iv = zeros(1, size(cf_used,2));
%size(iv)
%iv(:) = idx(1,:,1,1);
%cf_used_sorted = cf_used(:,iv,:,:);

%[cf_used_cf, ~, cf_used_n, cf_used_i] = ind2sub(size(cf_used), find(cf_used(:,27,:,:) ~= 1));

%[cf_used_cf, ~, cf_used_n, cf_used_i] = ind2sub(size(cf_used_change), find(cf_used_change > 0) );
%max(cf_used_cf)
%max(cf_used_n)
%max(cf_used_i)

%cf_used = sub2ind(size(cf), cf_used_cf, ones(length(cf_used_cf), 1), 1)
%size(cf_used_cf)
%size( ones(length(cf_used_cf), 1) )
% diff_xy = diff(xy(:,:,:), 1, 3);
% %size(xy)
% %size(diff_xy)
% 
% figure(909);
% plot( squeeze(diff_xy(1,:,:))' );
% 
% 
% figure(910);
% clf reset; % Reset figure.
% hold on;
% title('5-period moving average of difference in x-coordinate'); % Add title
% %for n = 1:pref.N
% n=1
% scatter(1:pref.iterations-1, squeeze(diff_xy(n,1,:)), 'Filled', 'k');
% simple = tsmovavg(squeeze(diff_xy(n,1,:)), 's', 5, 1);
% plot(simple, 'Color', color(n,:));
% simple = tsmovavg(squeeze(diff_xy(n,1,:)), 's', 10, 1);
% plot(simple, 'Color', color(n+1,:));
% simple = tsmovavg(squeeze(diff_xy(n,1,:)), 's', 100, 1);
% plot(simple, 'Color', color(n+2,:));
% %end
% plot(repmat(mean(squeeze(diff_xy(n,1,:))), 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');
% hold off;
% 
% figure(911);
% clf reset; % Reset figure.
% hold on;
% title('5-period moving average of difference in y-coordinate'); % Add title
% %for n = 1:pref.N
% n=1
% scatter(1:pref.iterations-1, squeeze(diff_xy(n,2,:)), 'Filled', 'k');
% simple = tsmovavg(squeeze(diff_xy(n,2,:)), 's', 5, 1);
% plot(simple, 'Color', color(n,:));
% simple = tsmovavg(squeeze(diff_xy(n,2,:)), 's', 10, 1);
% plot(simple, 'Color', color(n+1,:));
% simple = tsmovavg(squeeze(diff_xy(n,2,:)), 's', 100, 1);
% plot(simple, 'Color', color(n+2,:));
% %end
% plot(repmat(mean(squeeze(diff_xy(n,2,:))), 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');
% hold off;
% 
% figure(911);
% autocorr( squeeze(diff_xy(4,1,101:end)) )
% 
% [theta, rho] = cart2pol( diff_xy(:,1,:), diff_xy(:,2,:) );
% 
% figure(912);
% clf reset; % Reset figure.
% hold on;
% title('5-period moving average of direction'); % Add title
% %for n = 1:pref.N
% n=1
% scatter(1:pref.iterations-1, squeeze(theta(n,1,:)), 'Filled', 'k');
% simple = tsmovavg(squeeze(theta(n,1,:)), 's', 5, 1);
% plot(simple, 'Color', color(n,:));
% simple = tsmovavg(squeeze(theta(n,1,:)), 's', 10, 1);
% plot(simple, 'Color', color(n+1,:));
% simple = tsmovavg(squeeze(theta(n,1,:)), 's', 100, 1);
% plot(simple, 'Color', color(n+2,:));
% %end
% plot(repmat(mean(squeeze(theta(n,1,:))), 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');
% hold off;
% 
% 
% figure(913);
% clf reset; % Reset figure.
% hold on;
% title('5-period moving average of distance'); % Add title
% %for n = 1:pref.N
% n=1
% scatter(1:pref.iterations-1, squeeze(rho(n,1,:)), 'Filled', 'k');
% simple = tsmovavg(squeeze(rho(n,1,:)), 's', 5, 1);
% plot(simple, 'Color', color(n,:));
% simple = tsmovavg(squeeze(rho(n,1,:)), 's', 10, 1);
% plot(simple, 'Color', color(n+1,:));
% simple = tsmovavg(squeeze(rho(n,1,:)), 's', 100, 1);
% plot(simple, 'Color', color(n+2,:));
% %end
% plot(repmat(mean(squeeze(rho(n,1,:))), 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');
% hold off;
% %ylim([0 0.12]);
% 
% 
% 

figure(912);
clf reset; % Reset figure.
hold on;
title('5-period moving average of x-coordinate'); % Add title
%for n = 1:pref.N
n=2
scatter(1:pref.iterations, squeeze(xy(n,1,:)), 'Filled', 'k');
simple = tsmovavg(squeeze(xy(n,1,:)), 's', 4, 1);
plot(simple, 'Color', color(n,:));
simple = tsmovavg(squeeze(xy(n,1,:)), 's', 16, 1);
plot(simple, 'Color', color(n+1,:));
simple = tsmovavg(squeeze(xy(n,1,:)), 's', 64, 1);
plot(simple, 'Color', color(n+2,:));
%end
plot(repmat(mean(squeeze(xy(n,1,:))), 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');
hold off;


figure(913);
clf reset; % Reset figure.
hold on;
title('5-period moving average of y-coordinate'); % Add title
%for n = 1:pref.N
n=2
scatter(1:pref.iterations, squeeze(xy(n,2,:)), 'Filled', 'k');
simple = tsmovavg(squeeze(xy(n,2,:)), 's', 4, 1);
plot(simple, 'Color', color(n,:));
simple = tsmovavg(squeeze(xy(n,2,:)), 's', 16, 1);
plot(simple, 'Color', color(n+1,:));
simple = tsmovavg(squeeze(xy(n,2,:)), 's', 64, 1);
plot(simple, 'Color', color(n+2,:));
%end
plot(repmat(mean(squeeze(xy(n,2,:))), 1, pref.iterations), 'Color', 'k', 'LineStyle', ':');
hold off;


state_used = squeeze(sum(cf2(:,1:13,:,:),3,'omitnan'));
state_used_total = sum(sum(state_used,1,'omitnan'),3);
figure(914);
stem(state_used_total);
title('Total number of states occurrences'); % Add title


state_used = squeeze(sum(cf2(:,1:13,:,:),3,'omitnan'));
state_used_total = sum(sum(state_used,1,'omitnan'),3);
figure(914);
stem(state_used_total);
title('Total number of states occurrences'); % Add title

