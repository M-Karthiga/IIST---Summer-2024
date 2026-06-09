function [t, u1, u2, u3, theta, phi] = simulateSpheroid(re, a, w1, w2, k2, theta0_deg, phi0_deg, t_end)
    Q = 0.5;
    k = [0, k2, 0];
    [P, R] = computePR(re);
    theta0 = deg2rad(theta0_deg);
    phi0 = deg2rad(phi0_deg);
    phi1_0 = phi0 * sin(theta0);   % φ₁ = φ * sinθ
    y0 = [theta0; phi1_0];

    params = struct('P', P, 'Q', Q, 'R', R, ...
                    'a', a, 'w1', w1, 'w2', w2, ...
                    'k', k);

    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t, y] = ode45(@(t, y) orientationODE(t, y, params), [0, t_end], y0, options);

    theta = y(:,1);
    phi1 = y(:,2);
    sin_theta = sin(theta);
    sin_theta(abs(sin_theta) < 1e-6) = 1e-6;  % avoid division by zero
    phi = phi1 ./ sin_theta;
    
    u1 = sin(theta) .* cos(phi);
    u2 = sin(theta) .* sin(phi);
    u3 = cos(theta);

    % --- Unit vector constraint ---
    norm_error = max(abs(u1.^2 + u2.^2 + u3.^2 - 1));
    if norm_error > 1e-3
        warning('Unit vector constraint violated! Max error: %e', norm_error);
    end
end

