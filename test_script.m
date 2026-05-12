clear; close all; clc;

%% Hypersonic Glide Vehicle State-Space Linearization Demo
% State:
%   x = [r; lambda; phi; v; gamma; psi]
%
%   r      = radial distance from Earth center, m
%   lambda = longitude, rad
%   phi    = latitude, rad
%   v      = speed, m/s
%   gamma  = flight-path angle, rad
%   psi    = heading angle, rad
%
% Control:
%   sigma = bank angle, rad
%
% Continuous model:
%   xdot = f(x, sigma)
%
% Linearized continuous model:
%   delta_xdot = A_k delta_x + B_k delta_sigma
%
% Discrete model:
%   delta_x[k+1] = F_k delta_x[k] + G_k delta_sigma[k]

%% Parameters

params.Re   = 6371.1e3;
params.mu   = 3.986004418e14;
params.rho0 = 1.225;
params.hs   = 6700;

params.S = 3.99;
params.m = 1400;

params.CL = 0.085;
params.CD = 0.015;


params_vec = [
    6371.1e3;
    3.986004418e14;
    1.225;
    6700;
    3.99;
    1400;
    0.085;
    0.015
];

dt = 1.0;
tFinal = 300;
tspan = 0:dt:tFinal;

%% Initial Condition

x0 = [
    params.Re + 80000;   % r, m
    deg2rad(10);          % longitude, rad
    deg2rad(10);          % latitude, rad
    5500;                % speed, m/s
    deg2rad(-5);         % flight-path angle, rad
    deg2rad(90)          % heading, rad
];

sigma_fun = @(t,x) deg2rad(20);

%% Simulate Nominal Trajectory

[t, Xtraj_ode] = ode45(@(t,x) hgv_dynamics(x, sigma_fun(t,x), params), tspan, x0);

%Xtraj = Xtraj_ode.';   % 6 x N

load('measurement_data.mat');
for timeStep = 1:size(hypersonicPositionData, 1)
    cartesianState = hypersonicPositionData(timeStep, :);
    Xtraj(:, timeStep) = return_relevant_state_vector(cartesianState);
end
t = (1:1:size(hypersonicPositionData, 1));
N = size(Xtraj,2);

%% Linearize and Discretize Along Trajectory

nx = 6;
nu = 1;

A_ct = zeros(nx,nx,N);
B_ct = zeros(nx,nu,N);
F_dt = zeros(nx,nx,N);
G_dt = zeros(nx,nu,N);
fval = zeros(nx,N);
sigma_traj = zeros(1,N);

[f_sym, A_sym, B_sym, X, U, P] = hgv_dynamics_symbolic_no_coriolis();

test(:,1) = Xtraj(:,1); 
for k = 1:N-1
    xk = Xtraj(:,k);
    %sigmak = sigma_fun(t(k), xk);
    %sigma_traj(k) = sigmak;
    sigmak = 0;
    % [A_ct(:,:,k), B_ct(:,:,k), fval(:,k)] = linearize_hgv_fd(xk, sigmak, params);
    % [F_dt(:,:,k), G_dt(:,:,k)] = c2d_exact(A_ct(:,:,k), B_ct(:,:,k), dt);

    A0 = double(subs(A_sym, [X; U; P], [xk; sigmak; params_vec]));
    B0 = double(subs(B_sym, [X; U; P], [xk; sigmak; params_vec]));

    test(:, k+1) = (diag(xk) + (A0))*ones(6,1);

    % F = F_dt(:,:,k); G = G_dt(:,:,k);
    % Fi = F \ eye(nx);
    % Mi = Fi'*Yu(:,:,k)*Fi;
    % Om = Mi*G*inv(Qi + G'*Mi*G);
end

%% Save Results

save('hgv_linearized_discrete_state_space.mat', ...
     't', 'Xtraj', 'sigma_traj', ...
     'A_ct', 'B_ct', 'F_dt', 'G_dt', 'fval', ...
     'params', 'dt');

disp('Saved hgv_linearized_discrete_state_space.mat');

%% Plot Results

fig = figure;
hold on;
plot(rad2deg(Xtraj(2,:)), rad2deg(Xtraj(3,:)), 'LineWidth', 1.5);
plot(rad2deg(test(2,:)), rad2deg(test(3,:)), 'LineWidth', 1.5);
xlabel('Longitude, deg');
ylabel('Latitude, deg');
title('Nominal HGV Ground Track');
legend('True Data', 'Predicted Data')
grid on;
saveas(fig, "ground_track.jpg");

fig = figure;
hold on;
plot(t, (Xtraj(1,:) - params.Re)/1000, 'LineWidth', 1.5);
plot(t, (test(1,:) - params.Re)/1000, 'LineWidth', 1.5);
xlabel('Time, s');
ylabel('Altitude, km');
title('Altitude vs Time');
legend('True Data', 'Predicted Data')
grid on;
saveas(fig, "alt.jpg");

fig = figure;
hold on;
plot(t, Xtraj(4,:), 'LineWidth', 1.5);
plot(t, test(4,:), 'LineWidth', 1.5);
xlabel('Time, s');
ylabel('Speed, m/s');
title('Speed vs Time');
legend('True Data', 'Predicted Data')
grid on;
saveas(fig, "speed.jpg");

fig = figure;
hold on;
plot(t, rad2deg(Xtraj(5,:)), 'LineWidth', 1.5);
plot(t, rad2deg(test(5,:)), 'LineWidth', 1.5);
xlabel('Time, s');
ylabel('Flight-Path Angle, deg');
title('Flight-Path Angle vs Time');
legend('True Data', 'Predicted Data')
grid on;
saveas(fig, "flight_path.jpg");

disp('Done.');

%% Local Functions

function dx = hgv_dynamics(x, sigma, params)
%HGV_DYNAMICS Generic 3-DOF spherical-Earth HGV dynamics.

    r     = x(1);
    phi   = x(3);
    v     = x(4);
    gamma = x(5);
    psi   = x(6);

    Re   = params.Re;
    mu   = params.mu;
    rho0 = params.rho0;
    hs   = params.hs;
    S    = params.S;
    m    = params.m;
    CL   = params.CL;
    CD   = params.CD;

    h = r - Re;

    rho = rho0 * exp(-h/hs);
    qbar = 0.5 * rho * v^2;

    g = mu / r^2;

    D = qbar * CD * S / m;
    L = qbar * CL * S / m;

    check_singularity(x);

    dx = zeros(6,1);

    dx(1) = v * sin(gamma);

    dx(2) = (v * cos(gamma) * sin(psi)) / (r * cos(phi));

    dx(3) = (v * cos(gamma) * cos(psi)) / r;

    dx(4) = -D - g * sin(gamma);

    dx(5) = (1/v) * ( ...
              L * cos(sigma) ...
            + ((v^2/r) - g) * cos(gamma));

    dx(6) = (1/v) * ( ...
              (L * sin(sigma)) / cos(gamma) ...
            + (v^2/r) * cos(gamma) * sin(psi) * tan(phi));
end

function [A, B, f0] = linearize_hgv_fd(x0, sigma0, params)
%LINEARIZE_HGV_FD Finite-difference continuous-time linearization.

    nx = length(x0);
    nu = 1;

    f0 = hgv_dynamics(x0, sigma0, params);

    A = zeros(nx,nx);
    B = zeros(nx,nu);

    eps_base = sqrt(eps);

    x_scale = max(abs(x0), 1);

    for i = 1:nx
        dx = zeros(nx,1);
        h = eps_base * x_scale(i);

        dx(i) = h;

        fp = hgv_dynamics(x0 + dx, sigma0, params);
        fm = hgv_dynamics(x0 - dx, sigma0, params);

        A(:,i) = (fp - fm) / (2*h);
    end

    h_sigma = eps_base * max(abs(sigma0), 1);

    fp = hgv_dynamics(x0, sigma0 + h_sigma, params);
    fm = hgv_dynamics(x0, sigma0 - h_sigma, params);

    B(:,1) = (fp - fm) / (2*h_sigma);
end

function [F, G] = c2d_exact(A, B, dt)
%C2D_EXACT Exact zero-order-hold discretization.
%
% Continuous:
%   delta_xdot = A delta_x + B delta_u
%
% Discrete:
%   delta_x[k+1] = F delta_x[k] + G delta_u[k]

    nx = size(A,1);
    nu = size(B,2);

    M = [
        A, B;
        zeros(nu,nx), zeros(nu,nu)
    ];

    Md = expm(M * dt);

    F = Md(1:nx, 1:nx);
    G = Md(1:nx, nx+1:nx+nu);
end

function check_singularity(x)
%CHECK_SINGULARITY Warns/errors near coordinate singularities.

    r     = x(1);
    phi   = x(3);
    v     = x(4);
    gamma = x(5);

    tol = 1e-9;

    if abs(r) < tol
        error('Singularity: r is near zero.');
    end

    if abs(v) < tol
        error('Singularity: v is near zero.');
    end

    if abs(cos(phi)) < tol
        error('Singularity: cos(phi) is near zero, latitude is near +/-90 deg.');
    end

    if abs(cos(gamma)) < tol
        error('Singularity: cos(gamma) is near zero, flight-path angle is near +/-90 deg.');
    end
end

function [f_sym, A_sym, B_sym, X, U, P] = hgv_dynamics_symbolic_no_coriolis()
%HGV_DYNAMICS_SYMBOLIC_NO_CORIOLIS
% Symbolic 3-DOF spherical-Earth HGV equations of motion without Coriolis.
%
% State:
%   X = [r; lambda; phi; v; gamma; psi]
%
% Control:
%   U = sigma
%
% Parameters:
%   P = [Re; mu; rho0; hs; S; m; CL; CD]

    syms r lambda_sym phi_sym v gamma_sym psi_sym real
    syms sigma_sym real
    syms Re mu rho0 hs S m CL CD real

    X = [
        r;
        lambda_sym;
        phi_sym;
        v;
        gamma_sym;
        psi_sym
    ];

    U = sigma_sym;

    P = [
        Re;
        mu;
        rho0;
        hs;
        S;
        m;
        CL;
        CD
    ];

    % Altitude
    h = r - Re;

    % Atmosphere
    rho = rho0 * exp(-h / hs);

    % Dynamic pressure
    qbar = sym(1/2) * rho * v^2;

    % Gravity
    g = mu / r^2;

    % Aerodynamic accelerations
    D = qbar * CD * S / m;
    L = qbar * CL * S / m;

    % Equations of motion
    f_sym = sym(zeros(6,1));

    f_sym(1) = v * sin(gamma_sym);

    f_sym(2) = (v * cos(gamma_sym) * sin(psi_sym)) / (r * cos(phi_sym));

    f_sym(3) = (v * cos(gamma_sym) * cos(psi_sym)) / r;

    f_sym(4) = -D - g * sin(gamma_sym);

    f_sym(5) = (1/v) * ( ...
                 L * cos(sigma_sym) ...
               + ((v^2 / r) - g) * cos(gamma_sym));

    f_sym(6) = (1/v) * ( ...
                 (L * sin(sigma_sym)) / cos(gamma_sym) ...
               + (v^2 / r) * cos(gamma_sym) * sin(psi_sym) * tan(phi_sym));

    % Continuous-time linearization
    A_sym = simplify(jacobian(f_sym, X));
    B_sym = simplify(jacobian(f_sym, U));
end


% Functions for converstion
function R = return_conversion_from_ECEF_to_NED(lat, lon)
R = [-sin(lat)*cos(lon), -sin(lat)*sin(lon),  cos(lat);
    -sin(lon),           cos(lon),          0;
    -cos(lat)*cos(lon), -cos(lat)*sin(lon), -sin(lat)];
end

% Functions to turn measurements into necessary values
function vehicleFixedPosition_m = calculate_vehicle_fixed_position(azimuth, elevation, range, currentLocation)
phi = azimuth;
th = elevation;
r = range;

currentLocation = currentLocation(:);
relCartesian = [r*cosd(th)*cosd(phi);...
    r*cosd(th)*sind(phi);...
    r*sind(th)];

vehicleFixedPosition_m = (currentLocation + relCartesian);
end

function [theta_lat, phi_lon, r_m] = calculate_llr(fixedPosition_m)
earth_radius_m = 6371.8*1000;
lla = ecef2lla(fixedPosition_m(:)');

theta_lat = lla(1);   % degrees
phi_lon = lla(2);   % degrees
r_m = (lla(3) + earth_radius_m);   % m
end

function [v_speed, gamma, psi] = calculate_velocity_values(fixedState)

r_m = fixedState(1:3);        % position (m)
v = fixedState(4:6);     % velocity (m/s)

r_norm = norm(r_m);
v_norm = norm(v);

gamma_rad = asin(dot(r_m, v) / (r_norm * v_norm)); % radians
gamma = rad2deg(gamma_rad);

[theta_lat, phi_lon, ~] = calculate_llr(r_m);
rotation_matrix = return_conversion_from_ECEF_to_NED(deg2rad(theta_lat), deg2rad(phi_lon));

v_ned = rotation_matrix * v(:);

v_north = v_ned(1);
v_east  = v_ned(2);

heading = atan2(v_east, v_north); % radians
psi = mod(rad2deg(heading), 360);

psi = deg2rad(psi);
gamma = deg2rad(gamma);

v_speed = v_norm;

end

function hypersonicStateVector = return_relevant_state_vector(cartesianStateVector)

currFixedPosition_m = cartesianStateVector(1:3);
[theta_lat, phi_lon, r_m] = calculate_llr(currFixedPosition_m);
[v_speed, gamma, psi] = calculate_velocity_values(cartesianStateVector);

hypersonicStateVector = [r_m; deg2rad(theta_lat); deg2rad(phi_lon); v_speed; gamma; psi];

end
