% Task two of four
% simulate the AR(1) process 

% Simulate AR(1) process for case 3
T = 1000; % Number of observations
c = 1;    % Intercept term
phi1 = 1; % Coefficient
y0 = 0;   % Initial value
series = simulateAR1(T, c, phi1, y0);


% Estimate rho, standard error, and test statistic
[rho_hat, s_error, test_stat] = estimate_case3(series);

% Display results
disp(['Estimated rho: ', num2str(rho_hat)]);
disp(['Standard error: ', num2str(s_error)]);
disp(['Test statistic (rho_hat - 1) / s.e.(rho_hat): ', num2str(test_stat)]);

%% SImulate it 10000 times 
% Simulation and estimation for 10,000 iterations
n_simulations = 10000; % Number of simulations
T = 1000; % Number of observations
c = 1; % Intercept term
phi1 = 1; % Coefficient
y0 = 0; % Initial value

% Preallocate results
rho_hat_vec = zeros(n_simulations, 1);
s_error_vec = zeros(n_simulations, 1);
test_stat_vec = zeros(n_simulations, 1);

% Loop through simulations
for i = 1:n_simulations
    % Simulate AR(1) process for case 3
    series = simulateAR1(T, c, phi1, y0);
    
    % Estimate rho, standard error, and test statistic
    [rho_hat_vec(i), s_error_vec(i), test_stat_vec(i)] = estimate_case3(series);
end

% Display summary statistics
disp(['Mean rho_hat: ', num2str(mean(rho_hat_vec))]);
disp(['Mean standard error: ', num2str(mean(s_error_vec))]);
disp(['Mean test statistic: ', num2str(mean(test_stat_vec))]);

%% %Estimate kernel density of (rho_hat - 1) / s.e.(rho_hat)
[estimated_density, x_values] = ksdensity(test_stat_vec);

% Construct sequence from -5 to 4.9 with increments of 0.1
seq = -5:0.1:4.9;

% Evaluate the standard normal density at these points
normal_density = normpdf(seq, 0, 1); % Mean = 0, Std = 1

% Plot both densities
figure;
hold on;
plot(x_values, estimated_density, 'LineWidth', 1.5, 'DisplayName', 'Kernel Density of (\rhô - 1) / s.e.(\rhô)');
plot(seq, normal_density, '--', 'LineWidth', 1.5, 'DisplayName', 'Standard Normal Density');
legend('show');
xlabel('(\rhô - 1) / s.e.(\rhô)');
ylabel('Density');
title('Kernel Density vs. Standard Normal Density');
hold off;

%% Trend stationary process
% Define the simulation function for trend-stationary AR(1) process
function series = simulateTrendStationaryAR1(T, c, phi1, trend_coef, y0)
    series = zeros(T, 1);
    series(1) = y0;
    for t = 2:T
        series(t) = c + phi1 * series(t-1) + trend_coef * t + randn;
    end
end

% Set parameters
T = 100; % Number of observations
phi1 = 0.85; % AR(1) coefficient
c = 1; % Intercept term
trend_coef = 0.1; % Trend coefficient
y0 = 1 / (1 - phi1); % Starting value
n_simulations = 10000; % Number of simulations

% Preallocate results for the trend-stationary alternative
rho_hat_vec_trend = zeros(n_simulations, 1);
s_error_vec_trend = zeros(n_simulations, 1);
test_stat_vec_trend = zeros(n_simulations, 1);

% Perform simulations for the trend-stationary alternative
for i = 1:n_simulations
    % Simulate trend-stationary AR(1) process
    series = simulateTrendStationaryAR1(T, c, phi1, trend_coef, y0);
    
    % Estimate rho, standard error, and test statistic using case 3
    [rho_hat_vec_trend(i), s_error_vec_trend(i), test_stat_vec_trend(i)] = estimate_case3(series);
end

% Kernel density estimation for the trend-stationary alternative
[estimated_density_trend, x_values_trend] = ksdensity(test_stat_vec_trend);

% Overlay the kernel density for the null hypothesis (assume null data exists)
[estimated_density_null, x_values_null] = ksdensity(test_stat_vec); % Null test statistics

% Plot the densities
figure;
hold on;
plot(x_values_null, estimated_density_null, '--', 'LineWidth', 1.5, 'DisplayName', 'Null Hypothesis (\phi = 1)');
plot(x_values_trend, estimated_density_trend, 'LineWidth', 1.5, 'DisplayName', 'Trend-Stationary Alternative');
legend('show');
xlabel('(\rhô - 1) / s.e.(\rhô)');
ylabel('Density');
title('Kernel Densities: Null Hypothesis vs. Trend-Stationary Alternative');
hold off;

