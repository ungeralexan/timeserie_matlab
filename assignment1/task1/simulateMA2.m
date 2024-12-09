function y = simulateMA2(mu, theta1, theta2, T)

    epsilon = randn(T + 2, 1);  % [epsilon_{-1}, epsilon_0, epsilon_1, ..., epsilon_T]

    % Step 2: Preallocate y
    y = zeros(T, 1);

    % Step 3: Compute the MA(2) process
    for t = 1:T
        y(t) = mu + theta1 * epsilon(t + 1) + theta2 * epsilon(t) + epsilon(t + 2);
    end
end


% why T + 2: this is done so that t = 1 has history for epsilon t -1 and
% epsilon t-2 at the first time point
% therefore epsilon (t+2) corresponds to the current et
% epsilon (t+1) corresponds to thte et-1
% and epsilon(t) corresponds to the et
