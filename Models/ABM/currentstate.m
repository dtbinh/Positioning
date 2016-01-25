function J = currentstate(ti, xy, distance)
%CURRENTSTATE
%   Calculates a 13-bit descriptor of the state of each firm's position.
%   
%   Jonas K. Sekamane. 
%   Version 0.01

    N = size(distance, 2);
    JD = NaN(N, 13);

    % POSITION 1-5: Firm distance to population center greater than 0.1, 
    % 0.25, 0.4, 0.6, or 1.2 standard deviations.
    %distance = sqrt(sum(( xyti-pop_mu ).^2,2)); % centoriod distance
    fundamental_bins = [0.1 0.25 0.4 0.6 1.2];
    
    % POSITION 6-11: Greater than 5-period, 10-period or 15-period moving 
    % average.
    ma = [4 16 64]; % Number of periods in moving average.
%   ma_l = max( [ti-ma; ones(1, length(ma))] ); % The lower index should at least be 1.
    ma_l = max( [ti+1-ma; ones(1, length(ma))] ); % The lower index should at least be 1.
    ma_xy = NaN(N, 2, length(ma));
    for m = 1:length(ma)
        % Calculate simple moving average for each period
%       ma_xy(:,:,m) = mean( diff( xy(:, :, ma_l(m):ti), 1, 3), 3 );
        ma_xy(:,:,m) = mean( xy(:, :, ma_l(m):ti), 3 );
    end
    % Dealing with rounding errors by adding the term 1e3*eps.
    ma_xy = ma_xy + 1e3*eps;
    
%    diff_l = max(ti-1, 1); % The lower index should be at least 1.
    for n = 1:N
        % Firm distance to population center greater than ... 
        JD(n, 1:5) = ( distance(n) > fundamental_bins );
        
        % Reshape MA to form [x5 y5 x10 y10 x15 y15].
%         diff_xy = diff( xy(n,:,diff_l:ti), 1, 3);
%         if(ti==1)
%             diff_xy = zeros(1,2);
%         end
%         JD(n, 6:11) = reshape( (repmat(diff_xy,length(ma),1) > squeeze(ma_xy(n,:,:))' ), 1, length(ma)*2);
        JD(n, 6:11) = reshape( (repmat(xy(n,:,ti),length(ma),1) > squeeze(ma_xy(n,:,:))' ), 1, length(ma)*2);
    end
    
    
    
    % POSITION 12-13: Always on and always off.
    JD(:,12) = 1;
    JD(:,13) = 0;
    
    % Output variables
    J = JD;
end