function y = simulateAR1(T, c, phi, y0)

    % Preallocate the output vector for efficiency
    y = zeros(T, 1);

    % Draw the innovations epsilon_t from a standard normal distribution
    epsilon = randn(T, 1);

    % Initialize the previous value with y0
    y_prev = y0;

    % Generate the AR(1) process recursively
    for t = 1:T
        y(t) = c + phi * y_prev + epsilon(t);
        y_prev = y(t);  % Update y_prev for the next iteration
    end
end



