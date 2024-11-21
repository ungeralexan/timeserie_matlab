rng(200)

specification = [
    0.85, 0.3, -0.2;  % a)
    1.0, -0.2, -0.2;  % b)
    0.4, 0.5, -0.5;];

% size(specification,1) rows of the specification vectors, how much
results = zeros(size(specification,1),100);

% this iterates over each stochastic process
for spec_idx = 1:size(specification,1);

    %phi then stores the coefficients for the current stochastic process

    phi = specification(spec_idx,:);

    F = [phi(1), phi(2), phi(3);
        1, 0, 0;
        0, 1, 0];

    F11_vector = zeros(1,100);

    F_current = F

    for j = 1:100
        if j > 1

            F_current = F_current * F

        end
        F11_vector(j) = F_current(1,1)

    end
    results(spec_idx,:) = F11_vector

end


figure 
plot(1:100, results(1,:), 'LineWidth',1.5)
title('First with the F matrix calc')

figure
plot(1:100, results(2,:), 'LineWidth',1.5)
title('Second with the F matrix calc')

figure
plot(1:100, results(3,:), 'LineWidth',1.5)
title('Third with the F matrix cal')


%% Now lets see whether the plot is same with the eigenvalue decomposition
specifications = [
    0.85, 0.3, -0.2;  % a)
    1.0, -0.2, -0.2;  % b)
    0.4, 0.5, -0.5;    % c)
    ];

results = zeros(size(specifications,1),100);

for spec_idx = 1:size(specifications,1)
    phi = specifications(spec_idx,:)

    F = [phi(1), phi(2), phi(3);
         1, 0, 0;
         0, 1, 0];


    F11_vector = zeros(1, 100);

    [T, Lambda] = eig(F);
    T_inv = inv(T);

    for j = 1:100
        Lambda_j = Lambda^j;
        F_j = T * Lambda_j * T_inv;
        F11_vector(j) = F_j(1, 1);  % Extract the (1,1) element
    end

    results(spec_idx,:) = F11_vector

end


figure
plot(1:100, results(3,:), 'LineWidth',1.5)

%%
% Create the simulation function
function poly_values = simulateevaluations(phi1, phi2, phi3, z_values)
    % Ensure that z_values is a vector
    if ~isvector(z_values)
        error('z_values must be a vector');
    end
    
    % Evaluate the polynomial for each value of z
    poly_values = 1 - phi1 .* z_values - phi2 .* z_values.^2 - phi3 .* z_values.^3;
end

% Create z_values and call the function
z_values = linspace(-2, 2, 500);
poly_values = simulateevaluations(0.85, 0.3, -0.2, z_values);

figure
plot(z_values, poly_values, 'LineWidth',1.5)
title('Polynomial Evaluation');
xlabel('z');
ylabel('Value');
hold on;
yline(0, '--','LineWidth',1.5);
hold off