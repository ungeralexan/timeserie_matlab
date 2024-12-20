%% Task 1 of assignmet 4

% Set simulation parameters
n_simulations = 10000;
T = 1000;

% Preallocate results for each case
cases = [1, 2, 4];
results_case = zeros(length(cases), n_simulations);

% Simulate and store results for each case
for c_idx = 1:length(cases)
    c = cases(c_idx);
    for i = 1:n_simulations
        if c == 1
            series = simulateAR1(T, 0, 1, 0); % Case 1
        elseif c == 2
            series = simulateAR1(T, 0, 0.5, 0); % Case 2
        elseif c == 4
            series = simulateAR1(T, 0, 0.9, 0); % Case 4
        end
        result = olsts(series);
        results_case(c_idx, i) = result(1); % Store only rho_hat
    end
end

% Plot Kernel Density Estimates
figure;
hold on;
for c_idx = 1:length(cases)
    rho_values = results_case(c_idx, :);
    [f, xi] = ksdensity(rho_values);
    plot(xi, f, 'LineWidth', 1.5, 'DisplayName', ['Case ' num2str(cases(c_idx))]);
    xline(mean(rho_values), '--', 'DisplayName', ['Mean Case ' num2str(cases(c_idx))]);
end
legend('show');
xlabel('rho estimate (\rhô)');
ylabel('Density');
title('Kernel Density Estimates for \rhô (Cases 1, 2, and 4)');
hold off;

% Print observations
for c_idx = 1:length(cases)
    disp(['Mean rho_hat for Case ', num2str(cases(c_idx)), ': ', num2str(mean(results_case(c_idx, :)))]);
end


% Got it The goal is to visualize the distribution of the 10,000 rhos
% estimates for each case of the dicky fuller test and compare the
% distributions using a kernel density plot. ​

%% Compute the test statistics for the three cases 

% Preallocate for test statistics
test_stat1 = zeros(length(cases), n_simulations); % T * (rho_hat - 1)
test_stat2 = zeros(length(cases), n_simulations); % (rho_hat - 1) / s.e.(rho_hat)

% Calculate test statistics for each case
for c_idx = 1:length(cases)
    for i = 1:n_simulations
        % Extract results from olsts for the current case and simulation
        rho_hat = results_case(c_idx, i);
        s_error = sqrt(results_case(c_idx, i));
        T_eff = T - 1; % Effective sample size
        
        % Compute test statistics
        test_stat1(c_idx, i) = T_eff * (rho_hat - 1);
        test_stat2(c_idx, i) = (rho_hat - 1) / s_error;
    end
end

% Plot Kernel Density for T(rho_hat - 1)
figure;
hold on;
colors = ['r', 'g', 'b']; % Define colors for the plots
for c_idx = 1:length(cases)
    [f, xi] = ksdensity(test_stat1(c_idx, :));
    plot(xi, f, 'LineWidth', 1.5, 'Color', colors(c_idx), 'DisplayName', ['Case \phi = ' num2str(cases(c_idx))]);
end
legend('show');
xlabel('T(\rhô - 1)');
ylabel('Density');
title('Kernel Density Estimates for T(\rhô - 1)');
hold off;

% Plot Kernel Density for (rho_hat - 1) / s.e.(rho_hat)
figure;
hold on;
for c_idx = 1:length(cases)
    [f, xi] = ksdensity(test_stat2(c_idx, :));
    plot(xi, f, 'LineWidth', 1.5, 'Color', colors(c_idx), 'DisplayName', ['Case \phi = ' num2str(cases(c_idx))]);
end
legend('show');
xlabel('(\rhô - 1) / s.e.(\rhô)');
ylabel('Density');
title('Kernel Density Estimates for (\rhô - 1) / s.e.(\rhô)');
hold off;


%% Task 3

% simulate my AR(1)
T = 100;
phi = 0.85;
c = 1;
y0 = 1/(1-phi); % stationary mean
n_simulations = 10000;

% Preallocate test statistics for stationary alternative
test_stat1_stationary = zeros(1, n_simulations);
test_stat2_stationary = zeros(1, n_simulations);

% Simulate and compute test statistics for stationary alternative
for i = 1:n_simulations
    % Simulate stationary AR(1) process
    series = simulateAR1(T, c, phi, y0);
    
    % Estimate rho using OLS
    output = olsts(series);
    rho_hat = output(1); % Extract rho estimate
    s_error = output(2); % Extract standard error of rho_hat
    
    % Compute test statistics
    T_eff = T - 1; % Effective sample size
    test_stat1_stationary(i) = T_eff * (rho_hat - 1);
    test_stat2_stationary(i) = (rho_hat - 1) / s_error;
end

% Plot kernel density for T(rho_hat - 1)
figure;
hold on;
[f_stationary1, xi_stationary1] = ksdensity(test_stat1_stationary);
plot(xi_stationary1, f_stationary1, 'LineWidth', 1.5, 'DisplayName', 'Stationary AR(1) T(\rhô - 1)');

% Add kernel density from null hypothesis for comparison (example for case 1)
[f_null1, xi_null1] = ksdensity(test_stat1(1, :));
plot(xi_null1, f_null1, 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Null Case 1 T(\rhô - 1)');

legend('show');
xlabel('T(\rhô - 1)');
ylabel('Density');
title('Kernel Density of T(\rhô - 1): Stationary vs. Null');
hold off;

% Plot kernel density for (rho_hat - 1) / s.e.(rho_hat)
figure;
hold on;
[f_stationary2, xi_stationary2] = ksdensity(test_stat2_stationary);
plot(xi_stationary2, f_stationary2, 'LineWidth', 1.5, 'DisplayName', 'Stationary AR(1) (\rhô - 1) / s.e.(\rhô)');

[f_null2, xi_null2] = ksdensity(test_stat2(1, :));
plot(xi_null2, f_null2, 'LineWidth', 1.5, 'LineStyle', '--', 'DisplayName', 'Null Case 1 (\rhô - 1) / s.e.(\rhô)');

legend('show');
xlabel('(\rhô - 1) / s.e.(\rhô)');
ylabel('Density');
title('Kernel Density of (\rhô - 1) / s.e.(\rhô): Stationary vs. Null');
hold off;

%%
% Simulated test statistic series (for demonstration)
% Replace with actual simulated values from your analysis
test_stat_null1 = randn(10000, 1) - 1; % Example null values for T(rho_hat - 1)
test_stat_null2 = randn(10000, 1);     % Example null values for (rho_hat - 1) / s.e.(rho_hat)
test_stat_alt1 = randn(10000, 1) - 0.5; % Example alternative values for T(rho_hat - 1)
test_stat_alt2 = randn(10000, 1) - 0.5; % Example alternative values for (rho_hat - 1) / s.e.(rho_hat)

% Significance level
alpha = 0.05;

% Calculate Type II error probabilities
[typeII_stat1, typeII_stat2] = calculateTypeIIError(test_stat_null1, test_stat_null2, ...
                                                    test_stat_alt1, test_stat_alt2, alpha);

% Display results
disp(['Type II error for T(rho_hat - 1): ', num2str(typeII_stat1)]);
disp(['Type II error for (rho_hat - 1) / s.e.(rho_hat): ', num2str(typeII_stat2)]);

% What is what:
%test_stat_null1: Simulated values of T(rho-1) under the null hypothesis of
%a unit root
%test_stat_null2: Simulated values of (rho -1)/standard error under the
%null hypothesis
%test_stat_alt1: Simulated values of T(rho-1) under the alternative
%hypothesis stationary AR(1)
% test_stat_alt2: Simulated values of (rho -1)/standard error under the
% alternative hypothesis 
%cd('path_to_your_folder')
%addpath('path_to_your_folder')