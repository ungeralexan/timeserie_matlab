function [typeII_error_stat1, typeII_error_stat2] = calculateTypeIIError(test_stat_null1, test_stat_null2, ...
                                                                        test_stat_alt1, test_stat_alt2, alpha)
    % Calculate critical values for the test statistics under the null hypothesis
    critical_value1 = prctile(test_stat_null1, alpha * 100); % Critical value for T(rho_hat - 1)
    critical_value2 = prctile(test_stat_null2, alpha * 100); % Critical value for (rho_hat - 1) / s.e.(rho_hat)
    
    % Calculate Type II error for T(rho_hat - 1)
    % Failing to reject the null hypothesis under the alternative
    typeII_error_stat1 = mean(test_stat_alt1 > critical_value1); 
    
    % Calculate Type II error for (rho_hat - 1) / s.e.(rho_hat)
    % Failing to reject the null hypothesis under the alternative
    typeII_error_stat2 = mean(test_stat_alt2 > critical_value2);
end
