%% ABM
% version 0.02
% Jonas K. Sekamane
%
% Exogenous number of firms.
%
% Inspired in part by: 
%   Lever and Sergenti (2011) and Arthur (2014, chapter 3).

clearvars;

%% 1. PREFERENCES
%pref.seed = rng('default');
%pref.seed = rng(17796749, 'twister');
pref.boundary = 10; % Number of standard deviations
pref.resolution = 50; % Length of the square (even number to include (0,0))
pref.N = 5; % Number of firms
pref.mu = .5; % Mean of subpopulation
pref.n_ratio = 2; % Relative size of subpopulation; n_l/n_r how much larger is the left subpopulation than the right subpopulation
pref.iterations = 200;
pref.M = 100; % Number of condition/forecast rules that each firm holds.
pref.a_a = 1-1/75; % Accuracy mememory paramenter.

% Decision rules
%pref.rules = repmat({'STICKER'},1,pref.N);
%pref.rules = repmat({'AGGREGATOR'},1,pref.N);
%pref.rules = repmat({'HUNTER'},1,pref.N);
pref.rules = repmat({'INDUCTOR'},1,pref.N);
%pref.rules = repmat({'PREDATOR'},1,pref.N);
%pref.rules = [{'PREDATOR'} {'STICKER'} {'STICKER'} {'STICKER'} {'STICKER'}];


%% 2. SETUP
%rng(pref.seed); % Seed such that the randomly generated results are repeatable
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
    % The left subpopulation share of total population
    weight = pref.n_ratio/(1+pref.n_ratio);
    % Population probability density function (PDF). Left subpopulation is 
    % left unchanged, while the right subpopulation will scale according to 
    % n_ratio. Deviding by 2^2 has no effect on results, but scales PDF so 
    % its comparable to the baseline distribution with mean (0,0) and 
    % std. dev. (1,1). The subpopulations halves the std dev., thus their 
    % variance will be one-fourth.
    F = (F_l + F_r*1/pref.n_ratio)/4;
    %F = F_l*weight + F_r*(1-weight); % Population pdf
    
    
    
    %%% 2.2 Population descriptive statistics
    % Mean and covariance -- Johnson (1987)
    
    pop_mu = mu_l*weight + mu_r*(1-weight); % Population mean along each dimension; [mu_x mu_y]
    pop_sigma = sigma_l*weight + sigma_r*(1-weight) + weight*(1-weight)*(mu_l-mu_r)*(mu_l-mu_r)'; % Population covariance matrix
    pop_sd_x = sqrt(pop_sigma(1,1)); % Population std. dev. along x-axis; sigma_x


    %%% 2.3 Firms
    % Randomy draw initial position of firms
    
    % start at (0,0) move uniformly up to 3 std. dev. in random direction
    [x_0, y_0] = pol2cart( rand(pref.N,1)*2*pi , rand(pref.N,1)*3 );
    xy_0 = [x_0 y_0];
    % bivariate normal distribution, mean 0, std. dev. 1, uncorrelated (rho=0)
    %xy_0 = normrnd(0,1,[2 pref.N])';
    % uniformly distributed on boundary
    %xy_0 = (rand(2, pref.N)*pref.boundary-pref.boundary/2)';
    
    % Initial condition/forecast rule
    J_l = 13; % Length of J / number of bits.

    % Conditions set to 1 and 0 each with probability 0.1, otherwise # / NaN.
    prob = [0.1 0.1 0.8];
    % Always 1 condition that matches everything, ie. all # / NaN. I use the 
    % first condition / row.
    cf_0_condition = NaN(pref.M, J_l, pref.N);
    % Creating probability bins.
    cum_prob = cumsum([0 prob]); 
    % M-1 remaning conditions. For each remaning condition draw uniform number [0,1] for each of the state in J.
    draw = rand(pref.M-1, J_l, pref.N); 
    % Give value 0, 1 or NaN depending on which bin the number between [0,1]
    % falls into. Numbers between 0.0-0.1 get value 0, between 0.1-0.2 get 1,
    % all other numbers get value NaN.
    for n=1:pref.N
        for c=2:pref.M
            cf_0_condition(c,:,n) = sum( bsxfun(@ge, draw(c-1,:,n)', cum_prob), 2)'-1;
        end
    end
    cf_0_condition(cf_0_condition==2) = NaN;

    % Intercept drawn uniformly from [-1.5, 1.5] std. dev. Each row of
    % intercepts have the form [B_x1, B_y1].
    cf_0_intercept = -1.5 + rand(pref.M, 2, pref.N) * 3;
    % Coeffecients along same dimension drawn uniformly from [-1.2, 1.2],
    % while coeffecients along opposing dimension drawn uniformly from 
    % [-0.2, 0.2]. Each row of betas have the form [B_xx B_xy B_yx B_yy].
    cf_0_coefficients(:,[1 4 3 2],:) = [-1.2 + rand(pref.M, 2, pref.N) * 2.4 ...
                                           -0.2 + rand(pref.M, 2, pref.N) * 0.4 ];

    % Covariance matrix is initially set such that the variance along same the 
    % dimension is 0.005 and the variance along the opposing dimension is 0. 
    % The form is [s2_xx s2_xy s2_yx s2_yy].
    cf_0_var = repmat([0.005 0 0 0.005], pref.M, 1, pref.N);

    % The initial accuracy of all condition/forecast is set at zero. 
    cf_0_acc = zeros(pref.M, 1, pref.N);

    cf_0 = [cf_0_condition cf_0_intercept cf_0_coefficients cf_0_var cf_0_acc];
    cf_i = cf_0;
    
    %%% 2.4 Prelocating / Preparing matrices for iteration
    
    % Creating a 3D matrix of firms. 
    % 1st dimension is firm, 2nd is (x,y) coordinates, and 3rd is iteration.
    xy                      = NaN(pref.N, 2, pref.iterations); 
    heading                 = NaN(pref.iterations, pref.N);  
    
    % Creating a 3D matrix of market and utility.
    % 1st dimension is x, 2nd is y, and 3rd is iteration.
    market = NaN(pref.resolution, pref.resolution, pref.iterations); % The cell/value is the closest firm
    utility = NaN(pref.resolution, pref.resolution, pref.iterations); % The cell/value is utility

    shares                  = NaN(pref.iterations, pref.N);
    rank                    = NaN(pref.iterations, pref.N);
    eccentricity            = NaN(pref.iterations, pref.N);
    mean_eccentricity       = NaN(pref.iterations, 1);
    ENP                     = NaN(pref.iterations, 1);
    %misery                  = NaN(pref.iterations, 1);
    %perimeter               = NaN(pref.iterations, pref.N);
    %perimeter_extrema       = NaN(pref.iterations, pref.N);
    centroid                = NaN(pref.N, 2, pref.iterations);
    centroid_distance       = NaN(pref.iterations, pref.N);
    cf                      = NaN([size(cf_i) pref.iterations]);
    J_states                = NaN(pref.N, 13, pref.iterations);
    cf_used                 = NaN(pref.N, size(cf_i, 2), pref.N-1, pref.iterations);

%% 3. EVOLUTION
h = waitbar(0,'Calculating evolution...');
for i = 1:pref.iterations
    waitbar(i/pref.iterations);
    if(i==1) 
        % For first iteration set initial firm position
        xy(:,:,i) = xy_0;
    else 
        
        % 13 bit descriptor of the current state of each firm.
        J = currentstate(i-1, xy, centroid_distance(i-1,:));
        
        active_cf_i = cell(pref.N,pref.N);
        
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
                        % if distance is less than 0.1 std. dev.
                        [dx_in, dy_in] = pol2cart(direction, 0.1);
                        if(rho>0.1)    
                            xy(n,:,i) = xy(n,:,i-1) + [dx_in dy_in];
                        else
                            xy(n,:,i) = xy(rank(i-1,1),:,i-1); % don't overshoot
                        end
                    else
                        xy(n,:,i) = xy(n,:,i-1);
                    end
                    
                
                %case 'EXPLORER'
                    % Survey several directions and move along the most profitable.

                
                case 'INDUCTOR'
                    % Locate to maximise market share but subject to the predicted movements of other firms.
 					% Predictions are generated by its most reliable hypothesis. 
					% Gradually discards poorly performing hypotheses and forms new hypotheses.
                    
                    [f_n_xy, f_n_rule, f_n_active] = forecast( n, xy(:,:,i-1), J, cf_i(:,:,n) );
                    
                    active_cf_i(:,n) = f_n_active;
                    cf_used(n, :, :, i) = cf_i(f_n_rule,:,n)';
                    
                    % Calculate market using forecast
%                    [~, f_n_utility_i] = marketshare4(f_n_xy, [X(:) Y(:)]);
%                    % New firm position based on the (predicted) most unsatisfied customer
%                    [~, f_n_position] = max(f_n_utility_i(:).*F(:)/sum(F(:)));
%                    xy_new = [X(f_n_position) Y(f_n_position)];
                    
                    % New firm position based on the (predicted) centroid
                    % of the delaunay triangle with largest market share.
                    xy_new = maxcov_delaunay(f_n_xy, [X(:) Y(:)], F);
                    
                    % Move directly to most unsatisfied customer position.
                    %xy(n,:,i) = xy_new;
                    
                    % Direction in direction of most unsatisfied customer
                    % position.
                    [direction, rho] = cart2pol(xy_new(:,1)-xy(n,1,i-1), xy_new(:,2)-xy(n,2,i-1));
                    % 0.1 std. dev. move in direction of most profitable competitor, 
                    % if distance is less than 0.1 std. dev.
                    [dx_in, dy_in] = pol2cart(direction, 0.1);
                    if(rho>0.1)    
                        xy(n,:,i) = xy(n,:,i-1) + [dx_in dy_in];
                    else
                        xy(n,:,i) = xy_new; % don't overshoot
                    end
                
            end
        end
        
        % Update accuracy
        cf_i = accuracy(xy(:,:,i), active_cf_i, cf_0, pref.a_a);
        
        % Save states
        J_states(:,:,i) = J;
        
    end
    
    %%% 3.2 Propeties
    
    % Save condition/dorecast rules
    cf(:,:,:,i) = cf_i;
    
    % Market and utility.
    [market_i, utility_i] = marketshare4(xy(:,:,i), [X(:) Y(:)]);
    market(:,:,i) = market_i;
    utility(:,:,i) = utility_i;
    
    % Market centroid coordinates
    F_firm_i  = NaN(pref.N,1);
    for firm=1:pref.N
        % Index for each customer of the firm.
        idx = find(market_i==firm);
        
        % Number of customers for each firm
        F_firm_i(firm) = sum( F(idx) );
        
        % Each xy-coordinat within the firm's market weighted with the probability density
        centroid(firm,:,i) = ([X(idx) Y(idx)]' * F(idx))' ./ F_firm_i(firm);
    end
    % Distance from each firm to the respective centroid of its market.
    centroid_distance(i,:) = sqrt(sum(( xy(:,:,i)-centroid(:,:,i) ).^2,2))';

    % Market shares
    shares(i,:) = F_firm_i/sum(F(:)); % Calculate market share for each firm
    [~,rank(i,:)] = sort(shares(i,:), 'descend'); % Ranking firms according to market share
    
    % Eccentricity 
    % Firm euclidean distance from center/mean of voter distribution).
    eccentricity(i,:) = sqrt(sum(( xy(:,:,i)-repmat(pop_mu,pref.N,1) ).^2,2)); % in std. dev.
    mean_eccentricity(i,:) = mean(eccentricity(i,:));
    
    % Effective number of firms (ENP) -- Lever and Sergenti (2011), Laakso and Taagepera (1979)
    ENP(i,:) = sum(F(:))^2 / sum(F_firm_i(:).^2);
    % Herfindahl-Hirschman Index (HHI) -- Eiselt and Marianov (2011)
    % HHI(i,:) = sum( ( F_firm/sum(F(:)) ).^2 ); % HHI = 1/ENP
    
    % Misery
    % Quadratic loss function / representativeness
    %misery(i,:) = sum(utility_i(:).*F(:))/sum(F(:));
  
    % Market perimeter
    %market_props_i = regionpropsext(market(:,:,i), 'Extrema', 'Perimeter');
%     perimeter(i,:) = cat(1,market_props_i.Perimeter)';
%     perimeter_extrema(i,:) = cat(1,market_props_i.ExtremaPerimeter)';
    
    
end
close(h);
