function dydt = orientationODE(t, y, p)
    theta = y(1);
    phi1  = y(2);


    P = p.P; Q = p.Q; R = p.R;
    a = p.a; w1 = p.w1; w2 = p.w2;
    k = p.k;

    sin_theta = sin(theta);
    if abs(sin_theta) < 1e-6
        sin_theta = 1e-6;
    end
    cos_theta = cos(theta);
    sin_2theta = sin(2 * theta);

    phi = phi1 / sin_theta;
    sin_2phi = sin(2 * phi);
    cos_2phi = cos(2 * phi);

    % --- ODEs from Eq. (12) in the paper ---
    dtheta = (P/2) * sin_2theta * sin_2phi * a * cos(w1 * t) ...
           + R * (cos_theta * cos(phi) * k(1) + cos_theta * sin(phi) * k(2) - sin_theta * k(3)) * cos(w2 * t);

    dphi1 = (P/2) * phi * cos_theta * sin_2theta * sin_2phi * a * cos(w1 * t) ...
          + sin_theta * (P * cos_2phi - Q) * a * cos(w1 * t) ...
          + R * (-sin(phi) * k(1) + cos(phi) * k(2) ...
          + phi * (cos_theta^2 * cos(phi) * k(1) + cos_theta^2 * sin(phi) * k(2) - sin_theta * cos_theta * k(3))) * cos(w2 * t);

    dydt = [dtheta; dphi1];
end


