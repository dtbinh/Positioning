
% mean of subpopulation
mu = 1.5;
% relative size of subpopulation
n_ration = 2; % n_l / n_r;

% Grid
x1 = -5:1:5;
x2 = -5:1:5;
[X1,X2] = meshgrid(x1,x2);

sd = 0.5; % standard deviation of each subpopulation

% Subpopulation right
mu_r = [mu 0]; % subpopulation mean only deviate on x-axis
sigma_r = [sd^2 0; 0 sd^2]; % subpopulation sd only deviate on x-axis. Uncorrelated bivariate distribution, pho=0;
F_r = mvnpdf([X1(:) X2(:)],mu_r,sigma_r);
F_r = reshape(F_r,length(x2),length(x1));

% Subpopulation left
mu_l = [-mu 0];
sigma_l = [sd^2 0; 0 sd^2];
F_l = mvnpdf([X1(:) X2(:)],mu_l,sigma_l);
F_l = reshape(F_l,length(x2),length(x1));

% Total population 
share = n_ration/(1+n_ration);
F = F_l*share + F_r*(1-share);
%F = F_l*(n_ration/(1+n_ration)) + F_r*(1/(1+n_ration));
%F = F_l*n_ration + F_r;

% 3D plot of density function
surf(x1,x2,F);%,'EdgeColor','none');
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-5 5 -5 5 0 .5])
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');

firm1 = [2.25683, 0.952177];
firm2 = [-1.98716, 3.75343];
firm3 = [-3.39219, -1.477];
firm4 = flip([-3.39219, -1.477]);
firm5 = [0, 0];
xy = [firm1; firm2; firm3; firm4];

% random firm locations
n = 5;
x0 = rand(1, n)*10-5; % Initial x-position of firm
y0 = rand(1, n)*10-5; % Initial y-position of firm
xy = [x0' y0'];

test

test2 = [-5 5; -2.5 5; -2.5 -5; -5 -5];
mvncdf([-5 -5],[-2.5 5],mu_l,sigma_l)

in = inpolygon(X1,X2, test2(:,1),test2(:,2));

sum(F_l(in))/sum(F_l(:))

% 2D plot of contours
contour(x1,x2,F,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
hold on;
scatter( xy(:,1)' , xy(:,2)', 'filled'); % Plot the firms with respective colors.
xlim([-6 6]); ylim([-6 6]);
xlabel('x'); ylabel('y');
%hold on;
%voronoi(xy(:,1)' , xy(:,2)');
plot(vx,vy);
plot(bbox(:,1),bbox(:,2));
scatter(xi, yi);
h = patch(test2(:,1),test2(:,2),3);
set(h,'EdgeColor','none','FaceAlpha',0.25);
hold off;

%h = voronoi(xy(:,1)' , xy(:,2)');


% [V,C] = voronoin(xy);
% 
for i = 1:length(C3)
     C3{i}
end

% for i = 1:length(C)
%     if all(C{i}~=1) % If at least one of the indices is 1,
%         % then it is an open region and we can't
%         % patch that.
%         h = patch(V(C{i},1),V(C{i},2),i); % use color i.
%         get(h)
%     end
% end

% voronoi(xy(:,1), xy(:,2));
% xlim([-5 5]); ylim([-5 5]);
% hold on;
% C5 = V(C{4},:);
% C5(end+1,:) = C5(1,:);
% plot(C5(:,1),C5(:,2),'.-');
% scatter(V(:,1),V(:,2));
% %xlim([-5 5]); ylim([-5 5]);
% hold off;




% population pdf
pop_F = reshape(F, flip(size(x1).*size(x2)));

% population mean
pop_mu = mu_l*share + mu_r*(1-share)

%sqrt(sigma_l^2 + sigma_r^2)

%pop_sigma = sigma_l*share^2 + sigma_r*(1-share)^2

%pop_sigma = (sigma_l + mu_l*mu_l')*share + (sigma_r + mu_r*mu_r')*(1-share) - pop_mu*pop_mu'

%median(F,2)

%sqrt(pop_sigma)
%C*(K1+m1*m1') + (C-1)*(K2+m2*m2') -m*m',

%sum(F(:))

%x_dim = x1.*mean(F,1);
%mu_x = x1*mean(F,1)'

%mu_y = x2*mean(F,2) % should be very close to zero

%((x1.^2).*mean(F,1))'-mu_x
%sum(((x1.^2).*mean(F,1))'-mu_x)
%sd_x = sqrt( sum( (x1-mu_x).^2 .*mean(F,1) ) )
%sd_x = sqrt( sum(((x1.^2).*mean(F,1))'-mu_x))

%sd_x = sqrt(1/length(x1) * sum((x_dim'-mu_x).^2) )

%dim(x1)
%test = (x1 .* mean(F,1))'

%size(F)

%mvncdf(F)
%mvncdf([-5 -5],[0 0],mu_l,sigma_l)

%test = ;
%size(test)
%integrand = @(x1, x2) reshape(F_l, size(x1));
%integral2(integrand,0,5,0,5)

%integrand2 = @(x, y) reshape(mvnpdf([(0,5), (0,5)],[0,0],[1,0;0,1]), size(x));
%integral2(integrand2,0,5,0,5)
%integrand2 = @(x, y) reshape(mvnpdf([x(:), y(:)],[0,0],[1,0;0,1]), size(x));
%integral2(integrand2,0,5,0,5)

%integrand2 = @(x, y) reshape(F, size(x));
%integral2(integrand2,0,5,0,5)

% test = @(x, y) size(x);
% 
% test(0,5,0,5)
% 
% sum(F(:))
% 
% numel(F)
% 
% 
% F_r = mvnpdf([X1(:) X2(:)],mu_r,sigma_r)
% 
% mvncdf([-5 5],mu_r,sigma_r)
% 
% 
% F_r(5:6,7)
% X1(5:6,7)
% X2(5:6,7)
% 71




%mvncdf(F)
%mvncdf([0 0],[1 1],mu,Sigma)
%mvncdf([-.1 -5],[.1 5],mu_r,sigma_r)
%mvncdf([-5 -.1],[5 .1],mu_r,sigma_r)
%mean(sum(F,1))
%image(F*10000)


% surf(x1,x2,F_l);
% hold on;
% surf(x1,x2,F_r);
% caxis([min(F_l(:))-.5*range(F_l(:)),max(F_l(:))]);
% axis([-3 3 -3 3 0 .4])
% xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
% hold off;

%mvncdf([0 0],[1 1],mu,Sigma);
%contour(x1,x2,F1,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
%hold on;
%contour(x1,x2,F2,[.0001 .001 .01 .05:.1:.95 .99 .999 .9999]);
%xlabel('x'); ylabel('y');
%hold off;





% random population
%sd_1 = 5;
%sd_2 = 5;
%y_mean1 = 0;
%y_mean2 = 0;

%x_mean2 = rand * 15;
%x_mean1 = 0 - x_mean2;

%subpop1 = 500000 + rand * 166667;
%subpop2 = 1000000 - subpop1;

