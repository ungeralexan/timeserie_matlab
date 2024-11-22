%% Clear the workspace


%% Task b)
%set a seed first
rng(888)



% Generate the first MA(2) 
y_a = simulateMA2(0, 0.7, 0.3, 100);



%Generate the second MA(2) process
y_b = simulateMA2(4, 5, -17, 100);


%% Plot the values for both parameter specifications

%plot the first figure
figure 
plot(1:100, y_a, '-b', 'LineWidth',1.5)
hold on;
plot(1:100, 0 * ones(100,1), 'r-','LineWidth',1.5);
hold off;
title('MA (2) Process with mu=0, phi1=0.7, phi2=0.3 and t = 100');
xlabel('time');
ylabel('values');


% for process b
figure 
plot(1:100, y_b, 'g-', 'LineWidth',1.5)
hold on;
plot(1:100, 4 * ones(100,1), 'r--','LineWidth',1.5);
hold off;
title('MA (2) Process with mu=4, phi1=5, phi2= -17 and t = 100');
ylabel('values')
xlabel('time')

%% Create the autocorrelations

%empirical autocorrelations for realizations a
[r_a, lags_a, bounds_a] = autocorr(y_a, 'NumLags', 1);
empirical_rho1_a = r_a(2);  % r_a(1) corresponds to lag 0

% theoretical autocorrelations
theoretical_rho1_a = (0.7+ 0.3*0.7) / (1 + 0.7^2 + 0.3^2);


% empricial autocorrelations for realizarions b
[r_b, lags_b, bounds_b] = autocorr(y_b, 'NumLags', 1);
empirical_rho1_b = r_b(2);  % r_a(1) corresponds to lag 0


theoretical_rho1_b = (5 + 5*(-17))/ (1+ 5^5 + (-17)^2);

fprintf('The empirical correlation for the realization a is: %f\n', empirical_rho1_a)
fprintf('The theoretical correlation for the realization a is: %f\n', theoretical_rho1_a)
fprintf('The empirical correlation for the realization b is: %f\n', empirical_rho1_b)
fprintf('The theoretical correlation for the realization b is: %f\n', theoretical_rho1_b)

difference_a = empirical_rho1_a - theoretical_rho1_a;
difference_b = empirical_rho1_b - theoretical_rho1_b;

fprintf('The difference between empirical and theoretcial autocorrelation for process a is: %f\n',difference_a)
fprintf('The difference between empirical and theoretcial autocorrelation for process b is: %f\n',difference_b)

% big difference in b small in a

%% Task 5: Compute Ensemble Means and Variances for 30 Realizations

rng(150)
%Define the number of realizations
num_realizations = 30;

% Initialize matrices to store realizations
% Each column represents a realization
y_a_matrix = zeros(100, num_realizations); % 100 Rows and and thirty and variables


% Simulate 30 realizations for Case a)
for i = 1:num_realizations
    y_a_matrix(:, i) = simulateMA2(0, 0.7, 0.4, 100);
end

ensemble_mean_a = mean(y_a_matrix, 2); % Mean across columns

ensemble_variance_a = var(y_a_matrix, 0, 2)

% what I calculate here is the mean along each row which corresponds to
% each time point
%ensemble mean is calculated by taking the mean of each row across the 30 realizations.


figure
plot(1:100, ensemble_mean_a, 'r-', 'LineWidth',1.5);
hold on;
plot(1:100, ensemble_variance_a, 'k-','LineWidth',1.5)
hold off;
title('Plot the ensemble mean and variance over time for the first paramteres')
xlabel('time')
ylabel('values')
legend('Mean (Red)', 'Variance (Black)', 'Location', 'best') % Add legend with custom labels



% now for the second parameter specification
y_b_matrix = zeros(100, num_realizations); % 100 Rows and and thirty and variables


% Simulate 30 realizations for Case a)
for i = 1:num_realizations
    y_b_matrix(:, i) = simulateMA2(4, 5, -17, 100);
end

ensemble_mean_b = mean(y_b_matrix, 2); % Mean across columns

ensemble_variance_b = var(y_b_matrix, 0, 2)


% now plot the second process
figure
plot(1:100, ensemble_mean_b, 'g-', 'LineWidth',1.5);
hold on;
plot(1:100, ensemble_variance_b, 'b-','LineWidth',1.5)
hold off;
title('Plot the ensemble mean and variance over time for the second paramteres')
xlabel('time')
ylabel('values')
legend('Mean (Green)', 'Variance (Blue)', 'Location', 'best') % Add legend with custom labels

%% Compute autocorrelations by using the thirty correlations

% Define the number of lags
num_lags = 99;

% Initialize vectors to store autocorrelations
auto_corr_a = zeros(num_lags, 1);



% Compute autocorrelations for Case a
for t = 1:num_lags
    % Extract y_t and y_{t+1} across all realizations
    y_t = y_a_matrix(t, :)';      % 30 x 1 vector
    y_tp1 = y_a_matrix(t+1, :)';  % 30 x 1 vector
    
    % Compute the correlation
    auto_corr_a(t) = corr(y_t, y_tp1);
end

figure;
plot(1:99, auto_corr_a, 'b-', 'LineWidth', 2);
xlabel('Lag');
ylabel('Autocorrelation');
title('First-Order Autocorrelations Across Lags - Case a');
legend('Autocorrelation', 'Location', 'best', 'FontSize', 12, 'Box', 'off');

% check whether the autocoreelations are dependent on time t


% what happens : for each time point t it gathers the values y_t from 1 to
% 30 and for y_t+1 from 1 to 30 as well.
% Now for each t from 1 to 99 it computes the correlations between y_t and
% y_t+1 across the 30 realizations this gives you the correaltion between
% the time point t and t+1. 
% the reult is then stored in a vector called auto_corr_b(t) that ranges
% over 99 values. 
% So what I wanted to do was to calculate the correlations between two time
% values over realizations over time, meaning I take the values which range
% over 30 and take the correlation between the first 30 values at time
% period t and the next thirty values at time period t+1
% du gehst praktisch die reihen runter von eine t bis t+1 und berechnest
% die correlatioen zwischen den 30 observations für reihe 1 und den 30
% observations für reihe 2 was t+1 (across the realizations) was dir dann
% die autocorrelationen gibt







%% Task 7
% All three moments look time independet thats why they are stationary
% Stationarity doesnt depend on the paramter values, for finite order MA
% process are always stationary. 