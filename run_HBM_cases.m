function run_HBM_cases(re_list, alpha_list, omega1_list)
    % Run HBM for multiple combinations of re, alpha, omega1.
    % Loops over all combinations, runs HBM, plots, and saves figures.
    for i = 1:length(re_list)
        re = re_list(i); % Scalar
        for j = 1:length(alpha_list)
            alpha = alpha_list(j);
            for m = 1:length(omega1_list)
                omega1 = omega1_list(m);
                fprintf('Running case: re=%.2f, alpha=%.2f, omega1=%.2f\n', re, alpha, omega1);
                HBM_spheroid_dynamics_with_jeffery(re, alpha, omega1);
                prefix = ''; if re < 1, prefix = 'o'; end
                filename = sprintf('%sa=%.1f,w1=%.1f.png', prefix, alpha, omega1);
                set(gcf, 'PaperPositionMode', 'auto');
                set(gcf, 'Position', [100, 100, 1200, 800]);
                print(filename, '-dpng', '-r300');
            end
        end
    end
end

function HBM_spheroid_dynamics_with_jeffery(re, alpha, omega1)
    % Combines Jeffery initial guess with HBM solver
    % Generates phase space plots, 3D trajectory
    N = 400;
    omega2 = 1.0;
    k = [0; 0; 0]; % No torque for simple shear
    n_harmonics = 15;
    lambda = 1e2;
    Ncycles = 15;
    samples_per_cycle = 800;

    % Shape parameters
    P = (re.^2 - 1) ./ (2 * (re.^2 + 1));
    Q = 0.5;
    if re > 1
        b = log(re + sqrt(re.^2 - 1)) ./ (re .* sqrt(re.^2 - 1));
    else
        b = acos(re) ./ (re .* sqrt(1 - re.^2));
    end
    R = (3 * re.^4 .* b .* (2 * re.^2 - 1) - 1) ./ (16 * pi * (re.^2 - 1) .* (re.^2 + 1));

    % Time grid
    steady_shear = (omega1 == 0);
    if steady_shear
        omega_basis = alpha / (re + 1/re);
    else
        omega_basis = omega1;
    end
    T = 2 * pi / omega_basis;
    t = linspace(0, T, N);

    % Initial guess
    x0 = jeffery_init_and_use(re, alpha, omega_basis, n_harmonics, Ncycles, samples_per_cycle, steady_shear);

    % Solve HBM
    lb = -inf(size(x0));
    ub = inf(size(x0));
    opts = optimoptions('lsqnonlin', 'Display', 'iter', 'MaxIterations', 1000, ...
                       'MaxFunctionEvaluations', 2e5, 'FunctionTolerance', 1e-8, ...
                       'StepTolerance', 1e-8, 'DiffMinChange', 1e-8, 'DiffMaxChange', 0.1);
    coeffs = lsqnonlin(@(x) residual_HBM_lsq(x, t, alpha, omega_basis, omega2, k, P, Q, R, n_harmonics, lambda, steady_shear), ...
                       x0, lb, ub, opts);

    % Reconstruct
    [u1, u2, u3, du1_dt, du2_dt, du3_dt] = reconstruct_time_series(coeffs, t, omega_basis, n_harmonics);

    % Check unit norm
    norm_u = sqrt(u1.^2 + u2.^2 + u3.^2);
    fprintf('Max unit norm error for re=%.2f, alpha=%.2f, omega1=%.2f: %e\n', re, alpha, omega1, max(abs(norm_u - 1)));

    % Plot
    plot_results(t, u1, u2, u3, du1_dt, du2_dt, du3_dt, omega_basis);

    % Energy fraction
    eta = harmonic_energy_fraction(coeffs, n_harmonics, 5);
    fprintf('Harmonic energy fraction (first 5 harmonics) for re=%.2f, alpha=%.2f, omega1=%.2f: %s\n', ...
            re, alpha, omega1, mat2str(eta, 4));
end

function x0 = jeffery_init_and_use(re, a, omega1, n_harmonics, Ncycles, samples_per_cycle, steady_shear)
    Tcycle = 2 * pi / omega1;
    tspan = linspace(0, Ncycles * Tcycle, Ncycles * samples_per_cycle + 1);
    P = (re.^2 - 1) ./ (2 * (re.^2 + 1));
    Q = 0.5;
    u0 = [0.1; 0.4; sqrt(max(0, 1 - 0.1^2 - 0.4^2))];
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [tt, UU] = ode45(@(t, u) jeffery_rhs(t, u, P, Q, a, omega1, steady_shear), tspan, u0, opts);
    discard_cycles = 5;
    idx_keep = tt >= discard_cycles * Tcycle;
    t_keep = tt(idx_keep);
    U_keep = UU(idx_keep, :);
    T = Tcycle;
    num_coeffs = 2 * n_harmonics + 1;
    coeffs_all = zeros(3, num_coeffs);
    for j = 1:3
        uj = U_keep(:, j)';
        a0 = mean(uj);
        c = zeros(1, num_coeffs);
        c(1) = a0;
        for h = 1:n_harmonics
            cosvec = cos(h * omega1 * t_keep');
            sinvec = sin(h * omega1 * t_keep');
            ac = (2 / T) * trapz(t_keep, uj .* cosvec);
            bc = (2 / T) * trapz(t_keep, uj .* sinvec);
            c(2 * h) = ac;
            c(2 * h + 1) = bc;
        end
        coeffs_all(j, :) = c;
    end
    x0 = [coeffs_all(1, :)'; coeffs_all(2, :)'; coeffs_all(3, :)'];
    fprintf('Initial guess built: %d coefficients (per variable = %d)\n', length(x0), num_coeffs);
    fprintf('u1 a0 = %.4g, u2 a0 = %.4g, u3 a0 = %.4g\n', coeffs_all(1,1), coeffs_all(2,1), coeffs_all(3,1));
end

function du = jeffery_rhs(t, u, P, Q, a, omega1, steady_shear)
    u1 = u(1); u2 = u(2); u3 = u(3);
    if steady_shear
        c = a;
    else
        c = a * cos(omega1 * t);
    end
    du1 = P * u2 * (1 - 2*u1^2) * c + Q * u2 * c;
    du2 = P * u1 * (1 - 2*u2^2) * c - Q * u1 * c;
    du3 = -2 * P * u1 * u2 * u3 * c;
    du = [du1; du2; du3];
end

function Rvec = residual_HBM_lsq(x, t, a, omega1, omega2, k, P, Q, Rcoeff, n, lambda, steady_shear)
    N = length(t);
    num_coeffs = 2 * n + 1;
    U = zeros(3, N);
    dU = zeros(3, N);
    for j = 1:3
        coeffs = x((j-1)*num_coeffs+1 : j*num_coeffs);
        u_j = coeffs(1) * ones(1, N);
        du_j = zeros(1, N);
        for h = 1:n
            u_j = u_j + coeffs(2*h) * cos(h * omega1 * t) + coeffs(2*h + 1) * sin(h * omega1 * t);
            du_j = du_j - h * omega1 * coeffs(2*h) * sin(h * omega1 * t) + ...
                          h * omega1 * coeffs(2*h + 1) * cos(h * omega1 * t);
        end
        U(j, :) = u_j;
        dU(j, :) = du_j;
    end
    RHS = zeros(3, N);
    for i = 1:N
        u1 = U(1, i); u2 = U(2, i); u3 = U(3, i);
        k1 = k(1); k2 = k(2); k3 = k(3);
        b1 = u3 * (u3 * k1 - u1 * k3) - u2 * (u1 * k2 - u2 * k1);
        b2 = u1 * (u1 * k2 - u2 * k1) - u3 * (u2 * k3 - u3 * k2);
        b3 = u2 * (u2 * k3 - u3 * k2) - u1 * (u3 * k1 - u1 * k3);
        if steady_shear
            cos1 = 1.0;
        else
            cos1 = cos(omega1 * t(i));
        end
        cos2 = cos(omega2 * t(i));
        RHS(1, i) = P * u2 * (1 - 2*u1^2) * a * cos1 + Q * u2 * a * cos1 + Rcoeff * b1 * cos2;
        RHS(2, i) = P * u1 * (1 - 2*u2^2) * a * cos1 - Q * u1 * a * cos1 + Rcoeff * b2 * cos2;
        RHS(3, i) = -2 * P * u1 * u2 * u3 * a * cos1 + Rcoeff * b3 * cos2;
    end
    Res = dU - RHS;
    Rvec = [];
    for j = 1:3
        r = Res(j, :);
        c = zeros(1, num_coeffs);
        c(1) = mean(r);
        for h = 1:n
            c(2*h) = 2/N * sum(r .* cos(h * omega1 * t));
            c(2*h + 1) = 2/N * sum(r .* sin(h * omega1 * t));
        end
        Rvec = [Rvec; c(:)];
    end
    for i = 1:N
        norm_sq = sum(U(:, i).^2);
        Rvec = [Rvec; lambda * (norm_sq - 1)];
    end
end

function [u1, u2, u3, du1_dt, du2_dt, du3_dt] = reconstruct_time_series(coeffs, t, omega1, n_harmonics)
    n_coeffs = 2 * n_harmonics + 1;
    if length(coeffs) ~= 3 * n_coeffs
        error('Coefficient vector has incorrect size: expected %d, got %d', 3 * n_coeffs, length(coeffs));
    end
    a0 = coeffs(1:3);
    a = zeros(3, n_harmonics);
    b = zeros(3, n_harmonics);
    for j = 1:3
        base = (j - 1) * n_coeffs;
        a(j, :) = coeffs(base + 2 : base + 1 + n_harmonics);
        b(j, :) = coeffs(base + 2 + n_harmonics : base + n_coeffs);
    end
    u1 = a0(1) * ones(size(t));
    u2 = a0(2) * ones(size(t));
    u3 = a0(3) * ones(size(t));
    du1_dt = zeros(size(t));
    du2_dt = zeros(size(t));
    du3_dt = zeros(size(t));
    for h = 1:n_harmonics
        u1 = u1 + a(1, h) * cos(h * omega1 * t) + b(1, h) * sin(h * omega1 * t);
        u2 = u2 + a(2, h) * cos(h * omega1 * t) + b(2, h) * sin(h * omega1 * t);
        u3 = u3 + a(3, h) * cos(h * omega1 * t) + b(3, h) * sin(h * omega1 * t);
        du1_dt = du1_dt - h * omega1 * a(1, h) * sin(h * omega1 * t) + h * omega1 * b(1, h) * cos(h * omega1 * t);
        du2_dt = du2_dt - h * omega1 * a(2, h) * sin(h * omega1 * t) + h * omega1 * b(2, h) * cos(h * omega1 * t);
        du3_dt = du3_dt - h * omega1 * a(3, h) * sin(h * omega1 * t) + h * omega1 * b(3, h) * cos(h * omega1 * t);
    end
end

function eta = harmonic_energy_fraction(coeffs, n_harmonics, m_harmonics)
    n_coeffs = 2 * n_harmonics + 1;
    eta = zeros(3, 1);
    for j = 1:3
        base = (j - 1) * n_coeffs;
        a0 = coeffs(base + 1);
        a = coeffs(base + 2 : base + 1 + n_harmonics);
        b = coeffs(base + 2 + n_harmonics : base + n_coeffs);
        E_total = 0.5 * a0^2 + 0.5 * sum(a.^2 + b.^2);
        E_partial = 0.5 * a0^2 + 0.5 * sum(a(1:min(m_harmonics,n_harmonics)).^2 + b(1:min(m_harmonics,n_harmonics)).^2);
        eta(j) = E_partial / E_total;
    end
end

function plot_results(t, u1, u2, u3, du1_dt, du2_dt, du3_dt,omega1)
    % Plot time series, phase space, and 3D trajectory
    figure;
    subplot(2, 3, 1);
    plot(t, u1, 'b-', t, u2, 'r-', t, u3, 'g-');
    title('Time Series of Orientation Components');
    xlabel('Time'); ylabel('u_i'); legend('u_1', 'u_2', 'u_3');
    grid on;

    subplot(2, 3, 2);
    plot(u1, du1_dt, 'b-');
    title('Phase Space: u_1 vs du_1/dt');
    xlabel('u_1'); ylabel('du_1/dt'); grid on;

    subplot(2, 3, 3);
    plot(u2, du2_dt, 'r-');
    title('Phase Space: u_2 vs du_2/dt');
    xlabel('u_2'); ylabel('du_2/dt'); grid on;


    subplot(2, 3, 4);
    plot(u3, du3_dt, 'r-');
    title('Phase Space: u_3 vs du_3/dt');
    xlabel('u_3'); ylabel('du_3/dt'); grid on;

    subplot(2, 3, 5);
    plot3(u1, u2, u3, 'b-', 'LineWidth', 1.5);
    title('3D Orientation Trajectory');
    xlabel('u_1'); ylabel('u_2'); zlabel('u_3');
    grid on; axis equal;

    % Plot Poincaré sections at t = n * 2pi/omega1
    T = 2 * pi / omega1;
    indices = find(mod(t, T) < 1e-3);
    if isempty(indices)
        warning('No Poincaré section points found. Check time grid or omega1.');
        return;
    end

    subplot(2,3,6);
    scatter3(u1(indices), u2(indices), u3(indices), 10, 'b', 'filled');
    title('Poincaré Section at t = n * 2\pi/\omega_1');
    xlabel('u_1'); ylabel('u_2'); zlabel('u_3');
    grid on; axis equal;
end