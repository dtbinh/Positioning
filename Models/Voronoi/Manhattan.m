
pref.N = 5;
[x0, y0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*5 );
xy = [x0 y0];

pref.boundary = 10; % Number of standard deviations
pref.resolution = 1000; % Length of the square (even number to include (0,0))
b = pref.boundary/2;
[x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
[X,Y] = meshgrid(x,y);

[~, marketcb] = pdist2(xy, [X(:) Y(:)], 'cityblock', 'Smallest', 1);
marketcb = reshape(marketcb, size(X));

[~, marketec] = pdist2(xy, [X(:) Y(:)], 'euclidean', 'Smallest', 1);
marketec = reshape(marketec, size(X));

fig1 = figure(1);
imagesc(x, y, marketcb);
hold on;
scatter(xy(:,1), xy(:,2), 'filled', 'MarkerFaceColor', 'red');
hold off;
saveas(fig1,'fig/Voronoi_manhattan.pdf')


fig2 = figure(2);
imagesc(x, y, marketec);
hold on;
scatter(xy(:,1), xy(:,2), 'filled', 'MarkerFaceColor', 'red');
hold off;
saveas(fig2,'fig/Voronoi_euclidean.pdf')