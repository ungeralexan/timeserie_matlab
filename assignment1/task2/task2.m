%% Clear everything
clear


%% Call the function for different parameter specifications

rng(44)
y_a = simulateAR1(100,0,1,0);
y_b = simulateAR1(100,0.7,1,0);
y_c = simulateAR1(100,4.5,0.1,5);
y_d = simulateAR1(100,4.4,-0.1,4);
y_e = simulateAR1(100,4.5,0.9,45);





%% Plot the realizations of the process

% compte the expected values
E_y_c = 4.5 / (1 - 0.1);
E_y_d = 4.4 / (1 - (-0.1));
E_y_e = 4.5 / (1 - 0.9);
% not possible for the first two specifications as phi is not smaller then
% 1


%first plot
figure 
plot(1:100, y_a, 'r-','LineWidth',1.5)
title('AR (1) process with c = 0, φ = 1, y0 = 0')
xlabel('time')
ylabel('values')



%second plot
figure 
plot(1:100, y_b, 'b-','LineWidth',1.5)
title('AR (1) process with c = 0.7, φ = 1, y0 = 0')
xlabel('time')
ylabel('values')


%third plot
figure
plot(y_c, 'LineWidth', 1.2);
hold on;
yline(E_y_c, '--', 'LineWidth', 1, 'Color', 'red'); 
title('Series y_c');
xlabel('Time');
ylabel('y_c');
legend('Realization', 'Expected Value', 'Location', 'best');



% fourth plot
figure;
plot(y_d, 'k-', 'LineWidth', 1.2);
hold on;
yline(E_y_d, '--', 'LineWidth', 1, 'Color', 'blue');
title('Series y_d');
xlabel('Time');
ylabel('y_d');
legend('Realization', 'Expected Value', 'Location', 'best');




% fifth plot

figure;
plot(y_e, 'LineWidth', 1.2);
hold on;
yline(E_y_e, '--', 'LineWidth', 1, 'Color', 'red');
title('Series y_e');
xlabel('Time');
ylabel('y_e');
legend('Realization', 'Expected Value', 'Location', 'best');

%smaller phis meaning less persistent meaning they convert back to the mean
%faster then any big phis that are more persistant

%% How does the choice of φ affect the behavior of the AR(1) process?

%{
- If |φ| is smaller then one we have a stationary process which would yield a
mean reverting process with finite mean and variance
- If |φ| = 1 we have a non stationary random walk behaviour meaning (we
have non mean reverting)
- if |φ| > 1 we encounter an explosive process which is unstable
- Moreover positive signs of φ create smooth trajectories and negative
signs create oscillatory trajectories 
%}


%% Compute autocorrelations of first order


% Calculate theoretical autocorrelation for first order
theoretical_acf_y_c = 0.1;  % phi for y_c
theoretical_acf_y_d = -0.1; % phi for y_d
theoretical_acf_y_e = 0.9;  % phi for y_e


%calculate the empirical autocorrelations
empirical_acf_y_c = autocorr(y_c,'NumLags', 1);
empirical_acf_y_d = autocorr(y_d,'NumLags', 1);
empirical_acf_y_e = autocorr(y_e,'NumLags', 1);


fprintf('Autocorrelation of first order for y_c:\n');
fprintf('Theoretical: %.2f\n', theoretical_acf_y_c);
fprintf('Empirical: %.2f\n\n', empirical_acf_y_c(2));  % Use the second value since the first is lag 0

fprintf('Autocorrelation of first order for y_d:\n');
fprintf('Theoretical: %.2f\n', theoretical_acf_y_d);
fprintf('Empirical: %.2f\n\n', empirical_acf_y_d(2));

fprintf('Autocorrelation of first order for y_e:\n');
fprintf('Theoretical: %.2f\n', theoretical_acf_y_e);
fprintf('Empirical: %.2f\n\n', empirical_acf_y_e(2));

%% Analyze the stationary properties


% Parameters for the AR(1) process y_c
T = 100;              % Time steps
c = 4.5;              % Constant term
phi = 0.1;            % AR(1) coefficient
y0 = 5;               % Initial value
num_realizations = 30; % Number of realizations

% Initialize matrix to store the realizations
realizations = zeros(T, num_realizations);

rng(130)
% Generate 30 realizations of the AR(1) process
for i = 1:num_realizations
    realizations(:, i) = simulateAR1(T, c, phi, y0);
end

% Compute ensemble mean and variance for each time step
ensemble_mean = mean(realizations, 2);  % Mean across the 30 realizations at each time step
ensemble_variance = var(realizations, 0, 2);  % Variance across the 30 realizations at each time step

% Plot ensemble mean and variance
figure;

% Plot ensemble mean
subplot(2, 1, 1);
plot(1:T, ensemble_mean, 'LineWidth', 1.5);
title('Ensemble Mean of 30 Realizations of AR(1) Process y_c');
xlabel('Time');
ylabel('Mean');

% Plot ensemble variance
subplot(2, 1, 2);
plot(1:T, ensemble_variance, 'LineWidth', 1.5);
title('Ensemble Variance of 30 Realizations of AR(1) Process y_c');
xlabel('Time');
ylabel('Variance');

sgtitle('Ensemble Mean and Variance for 30 Realizations of AR(1) Process y_c');



% now lets do this for the process y_d


% Parameters for the AR(1) process y_d
T = 100;              % Time steps
c = 4.4;              % Constant term
phi = -0.1;            % AR(1) coefficient
y0 = 4;               % Initial value
num_realizations = 30; % Number of realizations

% Initialize matrix to store the realizations
realizations = zeros(T, num_realizations);

rng(131)
% Generate 30 realizations of the AR(1) process
for i = 1:num_realizations
    realizations(:, i) = simulateAR1(T, c, phi, y0);
end

% Compute ensemble mean and variance for each time step
ensemble_mean = mean(realizations, 2);  % Mean across the 30 realizations at each time step
ensemble_variance = var(realizations, 0, 2);  % Variance across the 30 realizations at each time step

% Plot ensemble mean and variance
figure;

% Plot ensemble mean
subplot(2, 1, 1);
plot(1:T, ensemble_mean,'r- ', 'LineWidth', 1.5);
title('Ensemble Mean of 30 Realizations of AR(1) Process y_d');
xlabel('Time');
ylabel('Mean');

% Plot ensemble variance
subplot(2, 1, 2);
plot(1:T, ensemble_variance,'r-', 'LineWidth', 1.5);
title('Ensemble Variance of 30 Realizations of AR(1) Process y_d');
xlabel('Time');
ylabel('Variance');

sgtitle('Ensemble Mean and Variance for 30 Realizations of AR(1) Process y_d');




% now lets do this for the process y_e


% Parameters for the AR(1) process y_e
T = 100;              % Time steps
c = 4.5;              % Constant term
phi = 0.9;            % AR(1) coefficient
y0 = 45;               % Initial value
num_realizations = 30; % Number of realizations

% Initialize matrix to store the realizations
realizations = zeros(T, num_realizations);

rng(132)
% Generate 30 realizations of the AR(1) process
for i = 1:num_realizations
    realizations(:, i) = simulateAR1(T, c, phi, y0);
end

% Compute ensemble mean and variance for each time step
ensemble_mean = mean(realizations, 2);  % Mean across the 30 realizations at each time step
ensemble_variance = var(realizations, 0, 2);  % Variance across the 30 realizations at each time step

% Plot ensemble mean and variance
figure;

% Plot ensemble mean
subplot(2, 1, 1);
plot(1:T, ensemble_mean,'g- ', 'LineWidth', 1.5);
title('Ensemble Mean of 30 Realizations of AR(1) Process y_e');
xlabel('Time');
ylabel('Mean');

% Plot ensemble variance
subplot(2, 1, 2);
plot(1:T, ensemble_variance,'g-', 'LineWidth', 1.5);
title('Ensemble Variance of 30 Realizations of AR(1) Process y_e');
xlabel('Time');
ylabel('Variance');

sgtitle('Ensemble Mean and Variance for 30 Realizations of AR(1) Process y_e');




%% This is the other way lets let it open for now
% Analyze properties of stationary AR process by using simulations

num_realizations = 30;

% preallocate arrays for storing realizations
realizations = zeros(100, num_realizations);


% Generate 30 realizations for the AR(3) process
for i = 1:num_realizations
    realizations(:,i) = simulateAR1(100,4.5,0.1,5);
end

ensembl_mean_c = mean(realizations,2);

ensemble_variance_c = var(realizations, 0, 2);



figure
plot(1:100, ensembl_mean_c, '-r', 'LineWidth',1.5);
hold on;
plot(1:100, ensemble_variance_c, 'b-','LineWidth',1.5);
hold off;

title('Plotted the ensemble mean and variances')
ylabel('values')
xlabel('time')


%{
So what did we do here : for each time step (row) from 1:100 we computed
the mean across all thirty realizations, which gave us the ensemble mean
for that specific time step
As result we got a series of 100 mean values, representing the ensemble
mean over time
%}



%% Compute autocorrelations by using the thirty correlations
% lets start with y_3
rng(150)
% Parameters
T = 100;               % Time steps
num_realizations = 30; % Number of realizations

% Initialize matrix to store the realizations
realizations = zeros(T, num_realizations);

% Generate 30 realizations of the AR(1) process
for i = 1:num_realizations
    realizations(:, i) = simulateAR1(T, 4.5, 0.1, 5);  % Using c = 4.5, phi = 0.1, y0 = 5
end

% Compute first-order autocorrelation for each realization
first_order_acfs = zeros(num_realizations, 1);
for i = 1:num_realizations
    acf = autocorr(realizations(:, i), 'NumLags', 1);
    first_order_acfs(i) = acf(2);  % First-order autocorrelation (second element)
end

% Plot the first-order autocorrelations for each realization
figure;
plot(1:num_realizations, first_order_acfs, 'o-', 'LineWidth', 1.5);
hold on;
yline(0.1, '--r', 'LineWidth', 1.2);  % Plot theoretical value as a reference
title('First-Order Autocorrelations for 30 Realizations of AR(1) Process y_c');
xlabel('Realization Index');
ylabel('First-Order Autocorrelation');
legend('Empirical Autocorrelation', 'Theoretical Autocorrelation (\phi = 0.1)', 'Location', 'best');
hold off;



% lets do this for y_d
rng(151)
% Parameters
T = 100;               % Time steps
num_realizations = 30; % Number of realizations

% Initialize matrix to store the realizations
realizations = zeros(T, num_realizations);

% Generate 30 realizations of the AR(1) process
for i = 1:num_realizations
    realizations(:, i) = simulateAR1(T, 4.4, -0.1, 4);  % Using c = 4.5, phi = 0.1, y0 = 5
end

% Compute first-order autocorrelation for each realization
first_order_acfs = zeros(num_realizations, 1);
for i = 1:num_realizations
    acf = autocorr(realizations(:, i), 'NumLags', 1);
    first_order_acfs(i) = acf(2);  % First-order autocorrelation (second element)
end

% Plot the first-order autocorrelations for each realization
figure;
plot(1:num_realizations, first_order_acfs, 'o-', 'LineWidth', 1.5);
hold on;
yline(-0.1, '--r', 'LineWidth', 1.2);  % Plot theoretical value as a reference
title('First-Order Autocorrelations for 30 Realizations of AR(1) Process y_d');
xlabel('Realization Index');
ylabel('First-Order Autocorrelation');
legend('Empirical Autocorrelation', 'Theoretical Autocorrelation (\phi = -0.1)', 'Location', 'best');
hold off;





% lets do this for y_e
rng(152)
% Parameters
T = 100;               % Time steps
num_realizations = 30; % Number of realizations

% Initialize matrix to store the realizations
realizations = zeros(T, num_realizations);

% Generate 30 realizations of the AR(1) process
for i = 1:num_realizations
    realizations(:, i) = simulateAR1(T, 4.5, 0.9, 45);  % Using c = 4.5, phi = 0.1, y0 = 5
end

% Compute first-order autocorrelation for each realization
first_order_acfs = zeros(num_realizations, 1);
for i = 1:num_realizations
    acf = autocorr(realizations(:, i), 'NumLags', 1);
    first_order_acfs(i) = acf(2);  % First-order autocorrelation (second element)
end

% Plot the first-order autocorrelations for each realization
figure;
plot(1:num_realizations, first_order_acfs, 'o-', 'LineWidth', 1.5);
hold on;
yline(0.9, '--r', 'LineWidth', 1.2);  % Plot theoretical value as a reference
title('First-Order Autocorrelations for 30 Realizations of AR(1) Process y_e');
xlabel('Realization Index');
ylabel('First-Order Autocorrelation');
legend('Empirical Autocorrelation', 'Theoretical Autocorrelation (\phi = 0.9)', 'Location', 'best');
hold off;


%% The other approach






% Define parameters
T = 100;                % Number of time steps
num_realizations = 30;  % Number of realizations

% Initialize matrix to store the 30 realizations
y_c_matrix = zeros(T, num_realizations);

% Generate 30 realizations of the AR(1) process
for i = 1:num_realizations
    y_c_matrix(:, i) = simulateAR1(T, 4.5, 0.1, 5);
end

% Initialize a vector to store first-order autocorrelations for each time step
first_order_autocorr = zeros(T - 1, 1);  % T - 1 because we use y_t and y_{t+1}

% Compute first-order autocorrelations for each time step
for t = 1:T-1
    % Extract y_t and y_{t+1} across all realizations (first-order lag)
    y_t = y_c_matrix(t, :)';      % Column vector for time t across realizations
    y_tp1 = y_c_matrix(t+1, :)';  % Column vector for time t+1 across realizations
    
    % Compute the first-order autocorrelation
    first_order_autocorr(t) = corr(y_t, y_tp1);
end

% Display the first-order autocorrelation for each time step
disp('First-order autocorrelations across time steps:');
disp(first_order_autocorr);




%Now lets plot ot
figure;
plot(1:99, first_order_autocorr, 'b-', 'LineWidth', 2);
xlabel('Lag');
ylabel('Autocorrelation');
title('First-Order Autocorrelations Across Lags - Case a');
legend('Autocorrelation', 'Location', 'best', 'FontSize', 12, 'Box', 'off');
