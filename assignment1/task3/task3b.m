%% Task 3b and c
% Lets redo the loops of task 3

%% Task 3b

% create a vector that contains all the psis and therefore weights of
% previous shocks F^j(1,1)


rng(399)
% all the psis for the three stochastic processes
specifications = [
    0.85, 0.3, -0.2;  % a)
    1.0, -0.2, -0.2;  % b)
    0.4, 0.5, -0.5;    % c)
    ]

%results vector with 100 columns with each psi for each process
results = zeros(size(specifications,1),100)

for spec_idx = 1:size(specifications,1)
    phi = specifications(spec_idx,:)

    F = [phi(1), phi(2), phi(3);
         1, 0, 0;
         0, 1, 0];


    F11_vector = zeros(1, 100);

    F_current = F
    for j = 1:100
        if j>1
            F_current = F_current * F
        end
        F11_vector(j) = F_current(1,1)
    end

    results(spec_idx, :) = F11_vector
end


figure 
hold on
for spec_idx = 1:size(specifications,1)
    plot(1:100, results(spec_idx,:), 'LineWidth',1.5)
end
title('Evolution of F^j(1,1) for Different AR(3) Specifications');
xlabel('j (Number of Steps)');
ylabel('F^j(1,1)');
legend('Specification a)', 'Specification b)', 'Specification c)', 'Location', 'Best');
grid on;
hold off;

% those 100 times refer to the 100 times steps in the future, the F(11) for
% every j then tracks the contribution of the initial state to the first
% variable after j steps. It tells you how the system changes in the long
% term

   
%% Now for one process

rng(200)

% this matrix contains rows where each row contains different psi
% coefficients for different stochastic processes
specifications = [
    0.85, 0.3, -0.2];

% the spec_idx is the index of the current row in specifications, it allows
% the programm to process each set of AR coefficients one by one


% final matrix hat 100 cols for 100 time perios and one row for the process
% which mean that we consider 100 psis
results = zeros(size(specifications,1),100);

% this iterates over each specification row, takes in into consideration
% each stochastic process
% it indexes the cureent row being processed here: the first row

for spec_idx = 1:size(specifications,1)



    phi = specifications(spec_idx,:)
    % the coeffciients of the current process are stored here
    % is like a small matrix that in each row entails the coefficient of
    % each stochastic process(we basicall extract the coefficients)

    F = [phi(1) , phi(2), phi(3);
        1 ,0, 0;
        0, 1, 0];
    % jsut plugged them into the F matrix

    F11_vector = zeros(1,100);
    % compute the F^j(11) matrix

    F_current = F
    for j = 1:100
        if j > 1
            % here for j 1 no multiplication happens yet
            F_current= F_current * F
        end
        F11_vector(j) = F_current(1,1)
    end

    results (spec_idx,:) = F11_vector

end


% now lets plot it 

figure 
hold on
plot(1:100, results, 'LineWidth',1.5)
title('Evolution of f^j(1,1) for different AR(3) specifications')
xlabel('Number of steps')
ylabel('F^j(1,1')
grid on
hold of


