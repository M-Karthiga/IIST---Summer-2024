re = 0.6;
a = 1;
w1 = 0.1;
w2 = 1.0;
k2 = 6.0;
h0 = 25;     %polar    
phi0 = 75;    %azimuthal standard conventions only   
t_end = 15000;   % total simulation time
transient_cutoff = 5000; 

% Run simulation 
[t1, u1w, u2w, u3w, theta, phi] = simulateSpheroid(re, a, w1, w2, k2, h0, phi0, t_end);

if t1(end) <= transient_cutoff
    error('Simulation ended at t = %.2f, before transient_cutoff = %.2f', t1(end), transient_cutoff);
end

% Remove transient part
idx = t1 > transient_cutoff;
t1 = t1(idx);
u1w = u1w(idx);
u2w = u2w(idx);
u3w = u3w(idx);
phi = phi(idx);
theta=theta(idx);

% Compute Derivative for Phase Space
t_uniform = linspace(transient_cutoff, t1(end), 10000)';  % uniform time grid
u1w_uniform = interp1(t1, u1w, t_uniform);              % interpolate u1
du1w_uniform = gradient(u1w_uniform, t_uniform(2) - t_uniform(1));       % compute du1/dt

u2w_uniform = interp1(t1, u2w, t_uniform);              % interpolate u2
du2w_uniform = gradient(u2w_uniform, t_uniform(2) - t_uniform(1));       % compute du2/dt

u3w_uniform = interp1(t1, u3w, t_uniform);              % interpolate u3
du3w_uniform = gradient(u3w_uniform, t_uniform(2) - t_uniform(1));       % compute du3/dt

% % Compute Second Derivatives
% d2u2w_uniform = gradient(du2w_uniform, t_uniform(2) - t_uniform(1));
% d2u3w_uniform = gradient(du3w_uniform, t_uniform(2) - t_uniform(1));
% 
% figure;
% plot3(u2w_uniform, du2w_uniform, d2u2w_uniform, 'g');
% xlabel('u_2'); ylabel('du_2/dt'); zlabel('d^2u_2/dt^2');
% title('3D Phase Space of u_2');
% grid on; axis tight;
% view([45 30]);
% 
% figure;
% plot3(u3w_uniform, du3w_uniform, d2u3w_uniform, 'g');
% xlabel('u_3'); ylabel('du_3/dt'); zlabel('d^2u_3/dt^2');
% title('3D Phase Space of u_3');
% grid on; axis tight;
% view([45 30]);


% === 3D Delay Embedding: u1(t), u1(t+τ), u1(t+2τ) ===
period = 2*pi / w2;           % dominant period
tau = period / 4;                             
dt = mean(diff(t1));                          
tau_index = round(tau / dt); % convert time delay to index

max_index = length(u1w) - 2*tau_index;
u1w_t0 = u1w(1:max_index);
u1w_t1 = u1w(tau_index+1 : tau_index+max_index);
u1w_t2 = u1w(2*tau_index+1 : 2*tau_index+max_index);

%Poincare sections
% --- Parameters
T = 2*pi / w2;                          % stroboscopic forcing period
t_strobe = transient_cutoff : T : t1(end);

% % --- Interpolate theta, phi at stroboscopic times
phi_strobe = interp1(t1, phi, t_strobe);       % azimuth (in code)
theta_strobe = interp1(t1, theta, t_strobe);   % polar (in code)


figure('Name', 'Dynamics Summary', 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

% === 1. 3D Trajectory Attractor ===
subplot(2, 3, 1);
plot3(u1w, u2w, u3w, 'b');
xlabel('u_1'); ylabel('u_2'); zlabel('u_3');
title('3D Attractor'); grid on; view([50 50]);

% === 2. Phase Space u1 ===
subplot(2, 3, 2);
plot(u1w_uniform, du1w_uniform, 'b');
xlabel('u_1'); ylabel('du_1/dt');
title('Phase Space: u_1'); xlim([-0.3 0.3]); ylim([-0.5 0.5]);
axis square; grid on;

% === 3. Phase Space u2 ===
subplot(2, 3, 3);
plot(u2w_uniform, du2w_uniform, 'b');
xlabel('u_2'); ylabel('du_2/dt');
title('Phase Space: u_2'); xlim([-0.5 0.5]); ylim([-1 1]);
axis square; grid on;

% === 4. Phase Space u3 ===
subplot(2, 3, 4);
plot(u3w_uniform, du3w_uniform, 'b');
xlabel('u_3'); ylabel('du_3/dt');
title('Phase Space: u_3'); xlim([0 1]); ylim([-0.3 0.3]);
axis square; grid on;

% === 5. Delay Embedding ===
subplot(2, 3, 5);
plot3(u1w_t0, u1w_t1, u1w_t2, 'm.');
xlabel('u_1(t)'); ylabel('u_1(t+\tau)'); zlabel('u_1(t+2\tau)');
title('3D Delay Embedding'); grid on;

% === 6. Poincaré Section ===
subplot(2, 3, 6);
scatter(phi_strobe, theta_strobe, 8, 'filled');
xlabel('\phi (deg)'); ylabel('\theta (deg)');
title('Poincaré Section'); axis square; grid on;

