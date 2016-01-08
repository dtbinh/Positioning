function O = burnin(M, estimate, stddev)
%BRUNIN   finds the 1st iteration where the summary variable is within 
%         +/- 1 std. dev. of the estimate.
%   
%   Jonas K. Sekamane. 
%   Version 0.01

    [repetitions, iterations, runs] = size(M);
    burnin = zeros(runs,repetitions);
    for run=1:runs
        for rep=1:repetitions
            serie = M(rep,:,run);
            lower = estimate(rep,:,run) - stddev(rep,:,run);
            upper = estimate(rep,:,run) + stddev(rep,:,run);
%           M(rep,:,run) = (lower<serie & serie<upper);
            burnin(run,rep) = find(lower<serie & serie<upper, 1); % Find 1st occurrence withn +/- 1 std. dev. of estimate.
            %mean(M(rep,burnin(rep,:,run):end,run)) % pct of observations after burn-in that are within +/- 1 std. dev.
        end
    end
    
    O = burnin;

end