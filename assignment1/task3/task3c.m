%% Task 3c Define AR coefficients for each specification
specifications = [
    0.85,  0.3, -0.2;  % Part a
    1.0,  -0.2,  0.2;  % Part b
    0.4,  -0.5,  0.5   % Part c
];

for s = 1:size(specifications, 1)
    phi1 = specifications(s, 1);
    phi2 = specifications(s, 2);
    phi3 = specifications(s, 3);
    
    F = [phi1, phi2, phi3; 1, 0, 0; 0, 1, 0];
    [T, Lambda] = eig(F);
    T_inv = inv(T);
    
    F_11_series_eig = zeros(1, 100);
    
    for j = 1:100
        Lambda_j = Lambda^j;
        F_j = T * Lambda_j * T_inv;
        F_11_series_eig(j) = F_j(1, 1);
    end
    
    figure;
    plot(1:100, F_11_series_eig, 'LineWidth', 1.5);
    title(['F_{11} series using Eigenvalue Decomposition for Specification ' char('a' + s - 1)]);
    xlabel('j');
    ylabel('F_j(1,1)');
    grid on;
end

%% For one process 
clear

rng(100)

% for only one specification of the process
specifications = [
    0.85,  0.3, -0.2  % Part a
];

% get the three phis 
% and aigain use s as the indice to specify about which specification of
% the process we are talking

% here I explicitly define that we start with the first process and so on
% (s can be seen as indice to loop across different specifications)
for s = 1:size(specifications, 1)
    phi1 = specifications(s, 1);
    phi2 = specifications(s, 2);
    phi3 = specifications(s, 3);
    
    %construct the F matrix
    F = [phi1, phi2, phi3;
        1, 0, 0;
        0, 1, 0];


    % the eigenvalue decomposition starts, T matrix contains the
    % eigenvectors and the lambda the eigenvalues of the F matrix
    [T, Lambda] = eig(F);
    T_inv = inv(T);
    
    % get the 100 eigenvalues for the 100 lags
    F_11_series_eig = zeros(1, 100);
    

    %now instead of multiplying F repeatidly, we raise the diagonal matrix
    %to power j
    for j = 1:100
        Lambda_j = Lambda^j;
        % reconstuct the F-j by simply using the new lambda matrix
        F_j = T * Lambda_j * T_inv;

        %again extrcat the F(11) element to construct the matrix
        F_11_series_eig(j) = F_j(1, 1);
    end
    
    figure;
    plot(1:100, F_11_series_eig, 'LineWidth', 1.5);
    title(['F_{11} series using Eigenvalue Decomposition for Specification ' char('a' + s - 1)]);
    xlabel('j');
    ylabel('F_j(1,1)');
    grid on;
end


%% now again doing it with the way I know
clear

rng(499)
% Define multiple specifications (processes)
specifications = [
    0.85,  0.3, -0.2;  % Specification 1
    0.5,  -0.2,  0.1;  % Specification 2
    -0.4,  0.6, -0.3   % Specification 3
];

% Once again choose an indice of the stochastic process I am using
num_specs = size(specifications, 1);
num_lags = 100;  % Number of lags to calculate
results_eig = zeros(num_specs, num_lags);

% Loop over each process specification
for s = 1:num_specs
    % Extract AR(3) coefficients for the current process
    phi1 = specifications(s, 1);
    phi2 = specifications(s, 2);
    phi3 = specifications(s, 3);
    
    % Construct the F matrix
    F = [phi1, phi2, phi3;
         1,    0,    0;
         0,    1,    0];
    
    % Compute eigenvalues and eigenvectors
    [T, Lambda] = eig(F);
    T_inv = inv(T);
    
    % Compute the F^{j}(1,1) series using eigenvalue decomposition
    F_11_series_eig = zeros(1, num_lags);
    for j = 1:num_lags
        Lambda_j = Lambda^j;
        F_j = T * Lambda_j * T_inv;
        F_11_series_eig(j) = F_j(1, 1);  % Extract the (1,1) element
    end
    
    % Store the results in the matrix
    results_eig(s, :) = F_11_series_eig;
end

% Display the results matrix (optional)
disp(results_eig);


