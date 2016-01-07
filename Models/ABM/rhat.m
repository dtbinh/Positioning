function O = rhat(M)
%RHAT   calculate the R-hat statistics of matrix with only the second half
%   
%   Jonas K. Sekamane. 
%   Version 0.03

    [m, n, runs] = size(M);
    array = [];
    for run=1:runs
        A = M(:,:,run); % Matrix with on iterations and repetitions for a particular run.
        
        % R-hat statistics / potential scale reduction factor (PSRF) -- Brooks and Gelman (1998) 
        
        % Within chain variance:
        W = mean( var(A, 0, 2) );
        %mean_j = mean(A,2);
        %var_j = 1/(n-1) * sum( (A-repmat(mean_j,1,n)).^2, 2 );
        %W = 1/m * sum(var_j);
        
        % Between chain variance
        B = n/(m-1) * sum(( mean(A,2) - mean(mean(A,1)) ).^2);
        %mean_ij = mean(mean_j);
        %B = n/(m-1) * sum( (mean_j-mean_ij).^2 );
        
        % Estimated varinace
        % (last term accounts for sampling variance of estimator).
        V_hat = (n-1)/n * W + B/n;% + B/(m*n) ;

        R_hat = sqrt(V_hat/W);
        array = [array R_hat];
    end
    O = array;

end