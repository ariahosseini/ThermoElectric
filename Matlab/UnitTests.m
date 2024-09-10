% Define test cases
test_cases = {
    struct('l', [1, 1], 'nk', [10, 10, 10], 'uo', 1e-3, 'm_frac', [1, 1, 1], 'v_frac', 1e-4, 'ko', [0, 0, 0], 'del_k', [0.1, 0.1, 0.1], 'n', 100, 'kind', 1),
    struct('l', [2, 2], 'nk', [20, 20, 20], 'uo', 2e-3, 'm_frac', [0.5, 1, 1.5], 'v_frac', 2e-4, 'ko', [0, 0, 0], 'del_k', [0.05, 0.05, 0.05], 'n', 50, 'kind', 2),
    struct('l', [1.5, 1.5], 'nk', [15, 15, 15], 'uo', 1e-2, 'm_frac', [1, 0.8, 1.2], 'v_frac', 1e-3, 'ko', [0, 0, 0], 'del_k', [0.1, 0.1, 0.1], 'n', 75, 'kind', 3)
};

% Initialize figure for plotting
figure;
hold on;
colors = lines(length(test_cases)); % Use different colors for each case

for i = 1:length(test_cases)
    % Extract parameters from the current test case
    params = test_cases{i};
    
    % Calculate scattering rates
    [mag_kpoint, E, tau, N] = inf_triangle(params.l, params.nk, params.uo, params.m_frac, params.v_frac, params.ko, params.del_k, params.n, params.kind);
    
    % Plot results
    plot(mag_kpoint, tau, 'DisplayName', sprintf('Test Case %d (kind=%d)', i, params.kind), 'Color', colors(i, :));
end

% Configure plot
xlabel('Magnitude of k-point (k)');
ylabel('Scattering Time (tau)');
title('Scattering Time vs. Magnitude of k-point for Different Test Cases');
legend('show');
grid on;
hold off;
