%% Load data
pop_mu = [-0.5000 0];

xy1 = squeeze(xy(1,:,1:500))';
iterations = length(xy1);

% Alternative
%xy1 = [-2.2280:(2.2280)/(500-1):0; 1.6269:-(1.6269)/(500-1):0]';

% Add some randomness, so distance is not always 0.1 std dev.
xy1 = xy1 + rand(iterations, 2)*0.05;

% Change
dxy1 = xy1(1:end-1,:) - xy1(2:end,:);

% Polar coordinates of change
[theta, rho] = cart2pol( dxy1(:,1), dxy1(:,2) );
v1 = [theta rho];


%% Tests


% TEST: multivaraite linear regression. Estimating this periods XY
% coordinates using last periods XY coordinates as dependent variable.
Y = xy1(2:end,:);
X = [ones(iterations-1,1) xy1(1:end-1,:)];
[beta, Sigma] = mvregress( X, Y );



% TEST: Descriptor
ti = 351;
xyti = xy1(ti,:);

% 13-bit descriptor of the state of the firm position.
J = NaN(1, 13);

% POSITION 1-5: Firm distance to population center greater than 0.2, 
% 0.6, 1, 1.6, or 2.4 standard deviations.
distance = sqrt(sum(( xyti-pop_mu ).^2,2)); % centoriod distance
fundamental_bins = [0.2 0.6 1 1.6 2.4];
J(1:5) = (distance>fundamental_bins);

% POSITION 6-11: Greater than 5-period, 10-period or 15-period moving 
% average.
ma = [5 10 15]; % Number of periods in moving average.
ma_l = max( [ti-ma; ones(1, length(ma))] ); % No lower index below 1.
ma_xy = NaN(length(ma), 2);
for m = 1:length(ma)
    % Calculate simple moving average for each period
    ma_xy(m,:) = mean( diff( xy1( ma_l(m):ti, : ), 1 ), 1 );
end;
% Reshape to form [x5 y5 x10 y10 x15 y15].
J(6:11) = reshape( (repmat(xyti,length(ma),1) > ma_xy), 1, length(ma)*2);

% POSITION 12-13: Always on and always off.
J(12) = 1;
J(13) = 0;



% TEST: Condition/forecast rule
M = 5;

% Initiate condition/forecast rule

% Conditions set to 1 and 0 each with probability 0.1, otherwise # / NaN.
prob = [0.1 0.1 0.8];
% Always 1 condition that matches everything, ie. all # / NaN. I use the 
% first condition / row.
cf_init_condition = NaN(M, length(J));
% Creating probability bins.
cum_prob = cumsum([0 prob]); 
% M-1 remaning conditions. For each remaning condition draw uniform number [0,1] for each of the state in J.
draw = rand(M-1, length(J)); 
% Give value 0, 1 or NaN depending on which bin the number between [0,1]
% falls into. Numbers between 0.0-0.1 get value 0, between 0.1-0.2 get 1,
% all other numbers get value NaN.
for c=2:M
    cf_init_condition(c,:) = sum( bsxfun(@ge, draw(c-1, :)', cum_prob), 2)'-1;
end
cf_init_condition(find(cf_init_condition==2)) = NaN;

% Intercept drawn uniformly from [-1.5, 1.5] std. dev. Each row of
% intercepts have the form [B_x1, B_y1].
cf_init_intercept = -1.5+rand(M,2)*3;
% Coeffecients along same dimension is drawn uniformly from [-1.2, 1.2],
% while coeffecients along opposing dimension is drawn uniformly from 
% [-0.2, 0.2]. Each row of betas have the form is [B_xx B_xy B_yx B_yy].
cf_init_coefficients(1:M, [1 3 2 4] ) = [-1.2+rand(M,2)*2.4 -.2+rand(M,2)*.4];

% Covariance matrix is initially set such that the variance along same the 
% dimension is 0.005 and the variance along the opposing dimension is 0. 
% The form is [s2_xx s2_xy s2_yx s2_yy].
cf_init_var = repmat([0.005 0 0 0.005], M, 1);

% The initial accuracy of all condition/forecast is set at zero. 
cf_init_acc = zeros(M, 1);

cf = [cf_init_condition cf_init_intercept cf_init_coefficients cf_init_var cf_init_acc];



% TEST: Find the conditions that satisfy the current state J.

% Subtracting the conditions from the state J gives 0 if a particular
% condition is statisfied, and 1 or -1 if it is not satisfied. Summing over
% all the particular conditions reveals if any of the conditions are 
% unfulfilled. By ignoring NaN when summing makes NaN a wildcard character, 
% that will match both 0s and 1s in state J.
count_unfulfilled = sum( abs( repmat(J, M, 1)-cf(:,1:13) ), 2, 'omitnan' );
% Index of conditions that satisfy the current state J.
c_idx = find( count_unfulfilled == 0);




%% Plots

% Draw evolution
for i = 1:500
    figure(2)
    clf reset; % Reset figure.
    n_xy = xy1(1:i,:);
    line(n_xy(:,1), n_xy(:,2), 'Color', [.8 .8 .8]);    
    xlim([-3 1]); ylim([-1.5 2.5]);
    title(sprintf('iteration %d',i)); % Add title
    hold on; % Place following scatter/plots on top of image.
    scatter( xy1(i,1) , xy1(i,2), 'filled'); % Plot the firms with respective colors.
    hold off; % Don't place anymore plots on top of figure.
    pause(.02);
end




figure(3)
iter = 500;
h = feather(xy1(1:iter,1), xy1(1:iter,2));
xlim([-5 iter+5]); ylim([-2 2]);
for l=1:length(h)-1
    h(l).XData(3:5) = [];
    h(l).YData(3:5) = [];
end
set(gca,'Position',[.03 .08 .96 .9])


figure(4)
iter = 100;
h = feather(dxy1(1:iter,1), dxy1(1:iter,2));
%xlim([-5 iter+5]); ylim([-.2 .2]);
for l=1:length(h)-1
    h(l).XData(3:5) = [];
    h(l).YData(3:5) = [];
end
set(gca,'Position',[.03 .08 .96 .9])

figure(5);
autocorr(xy1(:,1));
figure(6);
autocorr(xy1(:,2));