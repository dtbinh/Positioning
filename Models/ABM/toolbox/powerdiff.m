function O = powerdiff(est, std, n)
%POWERDIFF   calculates the power of the two-sample t-test where H_0: estimate is equal estimate from adjacent grid, H_A: different
%            Using the estimate and standard deviation of our summary variable.
%            Using a significane level of 0.05.
%            Where n is the number of post-burnin iterations.
%
%   Jonas K. Sekamane. 
%   Version 0.03

    
    p0 = [est std];
    power_diff = NaN(length(est), 1);
    for run=2:length(est)
        power_diff(run) = sampsizepwr('t2', p0(run,:), p0(run-1,1), [], n, 'Alpha', 0.05, 'Ratio', 1);
    end
    
    O = power_diff;
end