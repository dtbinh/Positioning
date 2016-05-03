function E = extrema(mu, n_ratio)
%% Determine "height" at all extrema point

%mu = rand * 1.5; % mu in [0, 1.5]
%n_ratio = 1 + rand; % n_ratio in [1, 2]

% Load extrapolation function
load('extrema_20160424_203825_ss10000.mat');


% There are 3 extrema when:
% n_ratio < - ( sqrt(4 mu^2 - 1) - 2mu ) / ( sqrt(4 mu^2 - 1) + 2mu ) e^(4 mu sqrt(4 mu^2 - 1) )

n_ratio_test = - ( sqrt(4*mu^2-1) - 2 * mu) / ( sqrt(4*mu^2-1) + 2 * mu) * exp( 4 * mu * sqrt(4*mu^2-1) );
%real(n_ratio_test)
cond_bimodal = (n_ratio < n_ratio_test);


% Find extrema for (mu, n_ratio) point

x_s = nan;
x_r = nan;
if(cond_bimodal)
    x_s = x_s_interpol(mu, n_ratio);
    x_r = x_r_interpol(mu, n_ratio);
end
x_l = x_l_interpol(mu, n_ratio);


%[mu n_ratio x_s x_r x_l]

E = [x_s x_r x_l];
end


% Finding density mass at extrema

%% POPULATION DISTRIBUTION
% pref.boundary = 10; % Number of standard deviations
% pref.resolution = 50; % Length of the square (even number to include (0,0))
% pref.mu = mu; % Mean of subpopulation
% pref.n_ratio = n_ratio; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
% sd = 0.5; % Standard deviation of each subpopulation
% b = pref.boundary/2;
% 
% [x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
% [X,Y] = meshgrid(x,y);
% 
% % Subpopulation right
% mu_r = [pref.mu 0]; % Subpopulation mean only deviate on x-axis
% sigma_r = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
% F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r); % Subpopulation pdf evaluated at each point in grid/square
% F_r = reshape(F_r,length(y),length(x)); % Formated as grid
% 
% % Subpopulation left
% mu_l = [-pref.mu 0]; % Subpopulation mean only deviate on x-axis
% sigma_l = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
% F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l); % Subpopulation pdf evaluated at each point in grid/square
% F_l = reshape(F_l,length(y),length(x)); % Formated as grid
% 
% % Total population 
% %weight = pref.n_ratio/(1+pref.n_ratio);
% F = (F_l + F_r*1/pref.n_ratio)/4;
% 
% 
% 
% 
% 
% F_x_s = nan;
% F_x_r = nan;
% if(cond_bimodal)
%     F_x_s = (mvnpdf([x_s 0],mu_l,sigma_l) + mvnpdf([x_s 0],mu_r,sigma_r)*1/pref.n_ratio)/4;
%     F_x_r = (mvnpdf([x_r 0],mu_l,sigma_l) + mvnpdf([x_r 0],mu_r,sigma_r)*1/pref.n_ratio)/4;
% end
% F_x_l = (mvnpdf([x_l 0],mu_l,sigma_l) + mvnpdf([x_l 0],mu_r,sigma_r)*1/pref.n_ratio)/4;
% 
% 
% %[F_x_s F_x_r F_x_l]
% 
% table(mu, n_ratio, x_s, x_r, x_l, F_x_s, F_x_r, F_x_l, max(F(:)), ...
%     'VariableNames', {'mu' 'n_ratio' 'x_s' 'x_r' 'x_l' 'F_x_s' 'F_x_r' 'F_x_l' 'F_max'})
%   
