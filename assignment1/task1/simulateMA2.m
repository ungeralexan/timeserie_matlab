function y = simulateMA2(mu, theta1, theta2, T)

    epsilon = randn(T + 2, 1);  % [epsilon_{-1}, epsilon_0, epsilon_1, ..., epsilon_T]

    % Step 2: Preallocate y
    y = zeros(T, 1);

    % Step 3: Compute the MA(2) process
    for t = 1:T
        y(t) = mu + theta1 * epsilon(t + 1) + theta2 * epsilon(t) + epsilon(t + 2);
    end
end

