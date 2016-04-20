% Two-step procedure to find the location that maximises market area

% Space
pref.boundary = 10; % Number of standard deviations
pref.resolution = 50; % Length of the square (even number to include (0,0))
b = pref.boundary/2;
[x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
[X,Y] = meshgrid(x,y);
%Population
pref.mu = 1.5; % Mean of subpopulation
pref.n_ratio = 2; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
sd = 0.5; % Standard deviation of each subpopulation
mu_r = [pref.mu 0]; % Subpopulation mean only deviate on x-axis
sigma_r = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r); % Subpopulation pdf evaluated at each point in grid/square
F_r = reshape(F_r,length(y),length(x)); % Formated as grid
mu_l = [-pref.mu 0]; % Subpopulation mean only deviate on x-axis
sigma_l = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l); % Subpopulation pdf evaluated at each point in grid/square
F_l = reshape(F_l,length(y),length(x)); % Formated as grid
weight = pref.n_ratio/(1+pref.n_ratio);
F = (F_l + F_r*1/pref.n_ratio)/4;

% Uniform
%F = ones(pref.resolution);

%pref.N = 12;
%[x0, y0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*3 );
%xy = [x0 y0];
xy = [-2 -2; -1 -3; 1 -0.5]; %[-2 -2; -1 -3; 1 -1];
n = size(xy,1);


%% Voronoi and Delaunay
%[vx,vy] = voronoi(xy(:,1)' , xy(:,2)');

[p,c,h,l, brel] = voronoib2(xy);
lvx = l(:,[1 3])';
lvy = l(:,[2 4])';

DT = delaunayTriangulation(xy(:,1), xy(:,2));

mar = 0.1;%0.0000000000001;


%% Step #1: Find potential
% Consider all (the easily calculable cases). 


% A) Intersection points between Voronoi edges and Delaunay edges. 
% There is little change in market area and shape when moving along 
% delaunay edge, but significant changes in area and shape when moving
% along Voronoi edge.

dvxy = reshape(DT.Points(DT.edges',:)', 4, [])';
cross = lineSegmentIntersect( l, dvxy );

crossxy = [cross.intMatrixX(cross.intAdjacencyMatrix) cross.intMatrixY(cross.intAdjacencyMatrix)];

% B) Intersection points of 3 or more Voronoi edges
[cc, ccr] = DT.circumcenter;

potentialxy = [cc; crossxy];

potentn = size(potentialxy,1);
for i=1:potentn
    [potenti_p, potenti_c, potenti_h, ~, ~] = voronoib2([potentialxy(i,:); xy]);
    potentpx = potenti_p(potenti_c{1},1);
    potentpy = potenti_p(potenti_c{1},2);
    potenthx{i} = potentpx(potenti_h{1});
    potenthy{i} = potentpy(potenti_h{1});
end


% C) Pairing locations.
% If a point is neighbour to a boundary: Take exsisting Voronoi region and 
% at the point split it parallel the neighbouring boundaries.

% Cut parallel to boundary lines.
potentn2 = 0;
potentn2 = size(brel,1);
for j=1:potentn2
    brelhori = (brel(j,2) > 2);
    cutline = [xy(brel(j,1),:); xy(brel(j,1),:)+[1-brelhori brelhori]];
    cutside = 1-(brel(j,2)-(1+2*brelhori)) + (1+2*brelhori);
    ployXY = p(c{brel(j,1)},:);
    %ployXY(h{brel(j,1)},:)
    polycut = cutpolygon(ployXY(h{brel(j,1)},:), cutline, cutside);
    potenthx{potentn+j} = polycut(:,1);
    potenthy{potentn+j} = polycut(:,2);
    potentialxy(potentn+j,:) = [xy(brel(j,1),:)+[1-brelhori brelhori].* mar];
end

% Cut parallel to voronoi lines
d2 = reshape(DT.edges, 1, [])';
d2vxy = [dvxy; dvxy(:,[3 4 1 2])];
d2diffx = d2vxy(:,3)-d2vxy(:,1);
d2diffy = d2vxy(:,4)-d2vxy(:,2);
d2diffinvnorm = 1 ./ sqrt(d2diffx.^2+d2diffy.^2);

cutlines = [d2vxy(:,1:2) d2vxy(:,1:2)+ [-d2diffy d2diffx]] ;
cutsides = 1-(d2diffy > 0)+1;

ppn = potentn+potentn2;
potentn3 = size(d2vxy, 1);
potentialxy(ppn+1:ppn+potentn3,:) = d2vxy(:,1:2)+[d2diffx d2diffy] .* repmat(d2diffinvnorm, 1, 2) .* mar;
for k=1:potentn3
    cutline = reshape(cutlines(k,[1 3 2 4]), 2, 2);
    ployXY = p(c{d2(k)},:);
    polycut = cutpolygon(ployXY(h{d2(k)},:), cutline, cutsides(k));
    potenthx{ppn+k} = polycut(:,1);
    potenthy{ppn+k} = polycut(:,2);
end


% -----
% Calculate potential
pn = potentn+potentn2+potentn3;
potential = NaN(pn,1);
for l=1:pn
   idx = inpolygon(X(:), Y(:), potenthx{l}, potenthy{l});
   potential(l) = sum(F(idx));
end

 
%% PLOTS

%i=10;
for i=1:ppn+potentn3
figure(21);
scatter(xy(:,1)', xy(:,2)', 'filled');
hold on;
patch(potenthx{i}, potenthy{i}, i, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
scatter(potentialxy(i,1), potentialxy(i,2), 'rd');
xlim([-6 6]); ylim([-6 6]);
hold off;
pause(0.5)
end

figure(25);
scatter(xy(:,1)', xy(:,2)', 'filled');
hold on;
for i=[3 5 6]%[3 5 6 13 16]
patch(potenthx{i}, potenthy{i}, i, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%scatter(potentialxy(i,1), potentialxy(i,2), 'rd');
end
xlim([-6 6]); ylim([-6 6]);
hold off;


figure(22);
scatter(xy(:,1)', xy(:,2)', 'filled');
hold on;
scatter(potentialxy(:,1), potentialxy(:,2), 'rd');
xlim([-6 6]); ylim([-6 6]);
%text(potentialxy(:,1)+rand(pn,1)-0.5, potentialxy(:,2)+rand(pn,1)-0.5, cellstr(num2str(potential)) );
hold off;

mask = sqrt(X.^2 + Y.^2) < 3; 
X_mask = X(mask==1);
Y_mask = Y(mask==1);

figure(23);
Fi = scatteredInterpolant(potentialxy(:,1),potentialxy(:,2),potential);
zi = Fi(X,Y);
%mesh(x,y,zi);
%Z = peaks(X,Y);
meshc(x,y,zi.*mask);
%plot3(potentialxy(:,1),potentialxy(:,2),potential,'.','markersize',12);
%xlim([-6 6]); ylim([-6 6]);
%grid on;

figure(24);
imagesc(x,y,zi.*mask);
set(gca,'ydir', 'normal');
xlim([-6 6]); ylim([-6 6]);


%tri = delaunay(xy(:,1)' , xy(:,2)');
%figure(3);
%triplot(tri, xy(:,1)' , xy(:,2)');
%xlim([-6 6]); ylim([-6 6]);


%[V,R] = voronoiDiagram(DT);
%DT.Points(DT.edges',:)
figure(4);
triplot(DT);
xlim([-6 6]); ylim([-6 6]);

%figure(1);
%plot(vx,vy, 'b');
%xlim([-6 6]); ylim([-6 6]);
%hold on;
%s = scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
%hold off;

figure(2);
scatter( xy(:,1)' , xy(:,2)', 'filled');
hold on;
%if(n>1) plot( all_split(:,[1 3])', all_split(:,[2 4])', 'k'); end
for i = 1:n
    X = p(c{i},1);
    Y = p(c{i},2);
    patch(X(h{i}), Y(h{i}), i, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end
xlim([-6 6]); ylim([-6 6]);
hold off;

[V1,B1] = voronoin(xy);

testcut = reshape(cutlines(1,[1 3 2 4]), 2, 2);
figure(5);
clf;
plot(lvx,lvy, 'b');
xlim([-6 6]); ylim([-6 6]);
hold on;
scatter( xy(:,1)' , xy(:,2)', 'filled', 'b');
%scatter( crossxy(:,1)' , crossxy(:,2)', 'r');
plot(testcut(:,1),testcut(:,2));
%patch(testcut(:,1), testcut(:,2), 1)
hold off;

%pref.seed = rng(97586720, 'twister');
%cutpolygon('demo')