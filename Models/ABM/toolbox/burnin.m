function O = burnin(M, varargin)
%BURNIN(M, ESTIMATE, STDDEV)   finds the 1st iteration where the summary variable is within 
%                              +/- 1 std. dev. of the estimate.
%
%BURNIN(M)   finds the 1st iteration where the summary variable is equal to
%            the last value.
%
%   Jonas K. Sekamane. 
%   Version 0.01

    switch nargin
        case 1
            % BURNIN(M)
            
            
            [repetitions, iterations, runs] = size(M);
            burnin = NaN(runs,repetitions);
            for run=1:runs
                for rep=1:repetitions
                    burnin(run,rep) = find(M(rep,:,run) == M(rep,end,run), 1); % Find 1st occurrence that equals last value
                    if(burnin(run,rep) == iterations)
                        warning(['Burnin only happened in at the last iteration! Repertition ' num2str(rep) ' in run ' num2str(run) ' may not have run long enough to determine burnin. Consider increasing the number of iterations']);
                    end
                end
            end

            O = burnin;
        
            
        case 3
            % BURNIN(M, ESTIMATE, STDDEV)
            
            estimate = varargin{1};
            stddev = varargin{2};
            
            [repetitions, iterations, runs] = size(M);
            burnin = NaN(runs,repetitions);
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

end