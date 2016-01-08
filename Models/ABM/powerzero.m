function O = powerzero(est, std, n)
%POWERZERO   calculates the power of the t-test where H_0: estimate is equal zero, H_A: different from zero
%            Using the estimate and standard deviation of our summary variable.
%            Using a significane level of 0.05.
%            Where n is the number of post-burnin iterations.
%
%   Jonas K. Sekamane. 
%   Version 0.03

    
    p0 = [est std];
    power_zero = zeros(length(est), 1);
    for run=1:length(est)
        power_zero(run) = sampsizepwr('t', p0(run,:), 0, [], n, 'Alpha', 0.05);
    end
    
    O = power_zero;
end