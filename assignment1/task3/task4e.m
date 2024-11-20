% Specifications for the three AR(3) processes
specifications = [
    0.85,  0.3, -0.2;  % Specification 1
    0.5,  -0.2,  0.1;  % Specification 2
    -0.4,  0.6, -0.3   % Specification 3
];

% Define the range for evaluating the polynomial (real z-values)
z_values = linspace(-2, 2, 500);  % Real z-values from -2 to 2

% Loop over each specification
for s = 1:size(specifications, 1)
    % Extract AR(3) coefficients
    phi1 = specifications(s, 1);
    phi2 = specifications(s, 2);
    phi3 = specifications(s, 3);
    
    % Call the function from Task 3.4
    [poly_values, poly_roots] = task4d(phi1, phi2, phi3, z_values);
    
    % Plot the polynomial evaluations
    figure;
    plot(z_values, poly_values, 'LineWidth', 1.5);
    hold on;
    
    % Highlight the roots on the plot
    scatter(real(poly_roots), zeros(size(poly_roots)), 50, 'r', 'filled', 'DisplayName', 'Roots');
    
    % Customize the plot
    xlabel('z');
    ylabel('Polynomial Value');
    title(['Polynomial for Specification ' num2str(s)]);
    legend('Polynomial', 'Roots');
    grid on;
end


%%
clear
% Without the loops and only for one

% Define AR(3) coefficients for one specification
phi1 = 0.85;
phi2 = 0.3;
phi3 = -0.2;

% Define the range for evaluating the polynomial (real z-values)
z_values = linspace(-2, 2, 500);  % Real z-values from -2 to 2

% Call the function from Task 3.4
[poly_values, poly_roots] = task4d(phi1, phi2, phi3, z_values);

% Plot the polynomial evaluations
figure;
plot(z_values, poly_values, 'LineWidth', 1.5);
hold on;


% Highlight the roots on the plot
scatter(real(poly_roots), zeros(size(poly_roots)), 50, 'r', 'filled', 'DisplayName', 'Roots');

% Customize the plot
xlabel('z');
ylabel('Polynomial Value');
title('Polynomial for Single Specification');
legend('Polynomial', 'Roots');
grid on;