function J = currentstate(ti, xy, distance)
%CURRENTSTATE
%   Calculates a 13-bit descriptor of the state of each firm's position.
%   
%   Jonas K. Sekamane. 
%   Version 0.01

    N = size(distance, 2);
    JD = NaN(N, 13);

    % POSITION 1-5: Firm distance to population center greater than 0.2, 
    % 0.6, 1, 1.6, or 2.4 standard deviations.
    %distance = sqrt(sum(( xyti-pop_mu ).^2,2)); % centoriod distance
    fundamental_bins = [0.2 0.6 1 1.6 2.4];
    
    % POSITION 6-11: Greater than 5-period, 10-period or 15-period moving 
    % average.
    ma = [5 10 15]; % Number of periods in moving average.
    ma_l = max( [ti+1-ma; ones(1, length(ma))] ); % The lower index should at least be 1.
    ma_xy = NaN(N, 2, length(ma));
    for m = 1:length(ma)
        % Calculate simple moving average for each period
        ma_xy(:,:,m) = mean( xy(:, :, ma_l(m):ti), 3 );
    end
    % Dealing with rounding errors by adding the term 1e3*eps.
    ma_xy = ma_xy + 1e3*eps;
    
    for n = 1:N
        % Firm distance to population center greater than ... 
        JD(n, 1:5) = ( distance(n) > fundamental_bins );
        
        % Reshape MA to form [x5 y5 x10 y10 x15 y15].
        JD(n, 6:11) = reshape( (repmat(xy(n,:,ti),length(ma),1) > squeeze(ma_xy(n,:,:))' ), 1, length(ma)*2);
    end
    
    
    
    % POSITION 12-13: Always on and always off.
    JD(:,12) = 1;
    JD(:,13) = 0;
    
    % Output variables
    J = JD;
end