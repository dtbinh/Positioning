%% ABM
% version 0.02
% Jonas K. Sekamane
%
% Exogenous number of firms (and thus decision rules).
%
% Inspired in part by: 
%   Lever and Sergenti (2011)

clearvars;

%% 1. PREFERENCES
pref.seed = rng('default');
pref.boundary = 10; % Number of standard deviations
pref.resolution = 70; % Length of the square (even number to include (0,0))
pref.N = 12; % Number of firms
pref.mu = 1.5; % Mean of subpopulation
pref.n_ratio = 2; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
pref.iterations = 200;

% Decision rules
%pref.rules = repmat({'STICKER'},1,pref.N);
%pref.rules = repmat({'AGGREGATOR'},1,pref.N);
pref.rules = repmat({'HUNTER'},1,pref.N);
%pref.rules = repmat({'PREDATOR'},1,pref.N);
%pref.rules = [{'PREDATOR'} {'STICKER'} {'STICKER'} {'STICKER'} {'STICKER'}];


%% 2. SETUP
rng(pref.seed); % Seed such that the randomly generated results are repeatable
unit = pref.resolution/pref.boundary; % 1 standard deviation
firms = 1:pref.N; % Gives each firm an ID/name
color = linspecer(pref.N); % Generates a color for each of the firms that is evenly space out in / distinguishable.
sd = 0.5; % Standard deviation of each subpopulation
b = pref.boundary/2;


    %%% 2.1 Population
    
    % Grid / Square
    [x,y] = deal( -b : pref.boundary/(pref.resolution-1) : b ); % Square with center at (0,0) and with length pref.boundary
    [X,Y] = meshgrid(x,y);

    % Subpopulation right
    mu_r = [pref.mu 0]; % Subpopulation mean only deviate on x-axis
    sigma_r = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
    F_r = mvnpdf([X(:) Y(:)],mu_r,sigma_r); % Subpopulation pdf evaluated at each point in grid/square
    F_r = reshape(F_r,length(y),length(x)); % Formated as grid

    % Subpopulation left
    mu_l = [-pref.mu 0]; % Subpopulation mean only deviate on x-axis
    sigma_l = [sd^2 0; 0 sd^2]; % Subpopulation std. dev. Uncorrelated bivariate distribution, rho=0;
    F_l = mvnpdf([X(:) Y(:)],mu_l,sigma_l); % Subpopulation pdf evaluated at each point in grid/square
    F_l = reshape(F_l,length(y),length(x)); % Formated as grid

    % Total population 
    weight = pref.n_ratio/(1+pref.n_ratio); % The left subpopulation share of total population
    F = F_l*weight + F_r*(1-weight); % Population pdf

    
    %%% 2.2 Population descriptive statistics
    % Mean and covariance -- Johnson (1987)
    
    pop_mu = mu_l*weight + mu_r*(1-weight); % Population mean along each dimension; [mu_x mu_y]
    pop_sigma = sigma_l*weight + sigma_r*(1-weight) + weight*(1-weight)*(mu_l-mu_r)*(mu_l-mu_r)'; % Population covariance matrix
    pop_sd_x = sqrt(pop_sigma(1,1)); % Population std. dev. along x-axis; sigma_x


    %%% 2.3 Firms
    % Randomy draw initial position of firms
    
    % start at (0,0) move uniformly up to 3 std. dev. in random direction
    %xy_0 = (cos(rand(2,pref.N)*2*pi) .* (rand(2,pref.N) * 3))';
    % bivariate normal distribution, mean 0, std. dev. 1, uncorrelated (rho=0)
    xy_0 = normrnd(0,1,[2 pref.N])';
    % uniformly distributed on boundary
    %xy_0 = (rand(2, pref.N)*pref.boundary-pref.boundary/2)';
    
   
    %%% 2.4 Prelocating / Preparing matrices for iteration
    
    % Creating a 3D matrix of firms. 
    % 1st dimension is firm, 2nd is (x,y) coordinates, and 3rd is iteration.
    xy                      = zeros(pref.N, 2, pref.iterations); 
    heading                 = zeros(pref.iterations, pref.N);  
    
    % Creating a 3D matrix of market and utility.
    % 1st dimension is x, 2nd is y, and 3rd is iteration.
    market = zeros(pref.resolution, pref.resolution, pref.iterations); % The cell/value is the closest firm
    utility = zeros(pref.resolution, pref.resolution, pref.iterations); % The cell/value is utility

    shares                  = zeros(pref.iterations, pref.N);
    rank                    = zeros(pref.iterations, pref.N);
    eccentricity            = zeros(pref.iterations, pref.N);
    mean_eccentricity       = zeros(pref.iterations, 1);
    ENP                     = zeros(pref.iterations, 1);
    misery                  = zeros(pref.iterations, 1);
    perimeter               = zeros(pref.iterations, pref.N);
    perimeter_extrema       = zeros(pref.iterations, pref.N);
    centroid                = zeros(pref.N, 2, pref.iterations);
    centroid_distance       = zeros(pref.iterations, pref.N);

%% 3. EVOLUTION
h = waitbar(0,'Calculating evolution...');
for i = 1:pref.iterations
    waitbar(i/pref.iterations);
    if(i==1) 
        % For first iteration set initial firm position
        xy(:,:,i) = xy_0;
    else 
        
        %%% 3.1 Decision rules
        % Firm position determined by firm behaviour / decsion rule
        
        for n = 1:pref.N
            switch pref.rules{n}
                case 'STICKER'
                    % Stay at current position
                    xy(n,:,i) = xy(n,:,i-1);

                    
                case 'AGGREGATOR'
                    % Move towards the center of its current market area.
                    % Maintain position is no customers
                    if(shares(i-1,n) > 0)
                        xy(n,:,i) = centroid(n,:,i-1);
                    else
                        xy(n,:,i) = xy(n,:,i-1);
                    end
                
                    
                case 'HUNTER'
                    % Continue same direction, if previous move proved fruitful,
                    % otherwise head the oppersite direction.
                     if(i==2)
                        % 0.1 std. dev. move in same direction, away from from (0,0)
                        [direction, ~] = cart2pol(xy(n,1,i-1), xy(n,2,i-1));
                        [dx_in, dy_in] = pol2cart(direction, 0.1);
                     else
                         [direction, ~] = cart2pol(xy(n,1,i-1)-xy(n,1,i-2), xy(n,2,i-1)-xy(n,2,i-2));
                         if(shares(i-1,n) > shares(i-2,n))
                             % 0.1 std. dev. move in same direction as previous move
                             [dx_in, dy_in] = pol2cart(direction, 0.1);
                         else
                             % 0.1 std. dev. move in oppersite direction (90 degree + random 180 degree)
                             [dx_in, dy_in] = pol2cart(direction + pi/2 + rand * pi, 0.1);
                         end
                     end
                     heading(i-1,n) = direction;
                     xy(n,:,i) = xy(n,:,i-1) + [dx_in dy_in];
                
                     
                case 'PREDATOR'
                    % Move towards its most profitable competitor.
                    if(rank(i-1,1) ~= n)
                        % Direction of most profitable competitor
                        [direction, rho] = cart2pol(xy(rank(i-1,1),1,i-1)-xy(n,1,i-1), xy(rank(i-1,1),2,i-1)-xy(n,2,i-1));
                        % 0.1 std. dev. move in direction of most profitable competitor, 
                        % if distance is less than 1 std. dev.
                        [dx_in, dy_in] = pol2cart(direction, 0.1);
                        if(rho>1)    
                            xy(n,:,i) = xy(n,:,i-1) + [dx_in dy_in];
                        else
                            xy(n,:,i) = xy(rank(i-1,1),:,i-1); % don't overshoot
                        end
                    else
                        xy(n,:,i) = xy(n,:,i-1);
                    end
                    
                
                %case 'EXPLORER'
                    % Survey several directions and move along the most profitable.

                
                %case 'INDUCTOR'
                    % Locate to maximise market share but subject to the predicted movements of other firms.
 					% Predictions are generated by its most reliable hypothesis. 
					% Gradually discards poorly performing hypotheses and forms new hypotheses.
                
            end
        end
        
    end
    
    %%% 3.2 Propeties
    
    % Market and utility.
    [market_i, utility_i] = marketshare3(xy(:,:,i), [X(:) Y(:)]);
    market(:,:,i) = market_i;
    utility(:,:,i) = utility_i;
    
    % Number of customers for each firm
    F_firm_i = arrayfun(@(firm) ...
            sum(F(find(market_i==firm))), ...
            firms);

    % Market shares
    shares(i,:) = F_firm_i/sum(F(:)); % Calculate market share for each firm
    [~,rank(i,:)] = sort(shares(i,:), 'descend'); % Ranking firms according to market share
    
    % Eccentricity 
    % Firm euclidean distance from center/mean of voter distribution).
    eccentricity(i,:) = pdist2(xy(:,:,i), pop_mu, 'euclidean'); % in std. dev.
    mean_eccentricity(i,:) = mean(eccentricity(i,:));
    
    % Effective number of firms (ENP) -- Lever and Sergenti (2011), Laakso and Taagepera (1979)
    ENP(i,:) = sum(F(:))^2 / sum(F_firm_i(:).^2);
    % Herfindahl-Hirschman Index (HHI) -- Eiselt and Marianov (2011)
    % HHI(i,:) = sum( ( F_firm/sum(F(:)) ).^2 ); % HHI = 1/ENP
    
    % Misery
    % Quadratic loss function / representativeness
    misery(i,:) = sum(utility_i(:).*F(:))/sum(F(:));
  
    % Market perimeter
    market_props_i = regionpropsext(market(:,:,i), 'Extrema', 'Perimeter');
%     perimeter(i,:) = cat(1,market_props_i.Perimeter)';
%     perimeter_extrema(i,:) = cat(1,market_props_i.ExtremaPerimeter)';
    
    % Market centroid coordinates
    for firm=1:pref.N
        % Index for each customer of the firm.
        idx = find(market_i==firm);
        % Each xy-coordinat within the firm's market weighted with the probability density
        centroid(firm,:,i) = ([X(idx) Y(idx)]' * F(idx))' ./ F_firm_i(firm);
    end
    % Distance from each firm to the respective centroid of its market.
    centroid_distance(i,:) = diag(pdist2(xy(:,:,i), centroid(:,:,i), 'Euclidean'))';
    
end
close(h);
