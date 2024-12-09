% Define the OLS estimation function for case 3
function [rho_hat, s_error, test_stat] = estimate_case3(series)
    T = length(series);
    y = series(2:T);          % Dependent variable
    X = [ones(T-1, 1), series(1:T-1)]; % Regressors: constant and lagged values
    beta = (X' * X) \ (X' * y); % OLS estimation
    rho_hat = beta(2);         % Extract rho estimate (slope)
    y_pred = X * beta;         % Predicted values
    residuals = y - y_pred;    % Residuals
    s_squared = sum(residuals.^2) / (T - 2); % Variance of residuals
    Var_rho = s_squared / sum((series(1:T-1) - mean(series(1:T-1))).^2); % Variance of rho_hat
    s_error = sqrt(Var_rho);   % Standard error of rho_hat
    test_stat = (rho_hat - 1) / s_error; % Test statistic
end