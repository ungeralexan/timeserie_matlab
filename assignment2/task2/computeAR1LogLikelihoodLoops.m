function logLikelihoodContributions = computeAR1LogLikelihoodLoops(params, y)
    % Extract parameters
    c = params(1);
    phi = params(2);
    sigma2 = params(3);

    % Number of log-likelihood contributions to compute
    T = length(y);
    logLikelihoodContributions = zeros(T-1, 1); % Preallocate for efficiency

    % Initialize y_prev
    y_prev = y(1);

    % Compute log-likelihood contributions iteratively
    for t = 2:T
        residual = y(t) - c - phi * y_prev; % Compute residual
        logLikelihoodContributions(t-1) = -0.5 * log(2 * pi * sigma2) - (residual^2) / (2 * sigma2);
        y_prev = y(t); % Update y_prev for the next iteration
    end
end

