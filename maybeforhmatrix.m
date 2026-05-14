function [z, H] = hgv_aer_speed_measurement(x, obs_ecef_m, obs_lat_rad, obs_lon_rad)
%HGV_AER_SPEED_MEASUREMENT
% Converts HGV state [r; lambda; phi; v; gamma; psi] to
% measurement [azimuth; elevation; range; speed] and returns linearized H.
%
% State:
%   x = [r; lambda; phi; v; gamma; psi]
%
% Inputs:
%   x           = 6x1 state vector
%                 r      = geocentric radius, m
%                 lambda = longitude, rad
%                 phi    = latitude, rad
%                 v      = speed, m/s
%                 gamma  = flight-path angle, rad
%                 psi    = heading angle, rad
%
%   obs_ecef_m  = 3x1 observer ECEF position, m
%   obs_lat_rad = observer latitude, rad
%   obs_lon_rad = observer longitude, rad
%
% Outputs:
%   z = [azimuth; elevation; range; speed]
%       azimuth   = rad, measured clockwise from north
%       elevation = rad
%       range     = m
%       speed     = m/s
%
%   H = dz/dx, 4x6 linearized measurement matrix

    % Extract state
    r      = x(1);
    lambda = x(2);
    phi    = x(3);
    v      = x(4);

    % Target ECEF position
    p_ecef = [
        r*cos(phi)*cos(lambda);
        r*cos(phi)*sin(lambda);
        r*sin(phi)
    ];

    % Line-of-sight vector in ECEF
    dp_ecef = p_ecef - obs_ecef_m(:);

    % ECEF to ENU rotation at observer
    A_enu = [
        -sin(obs_lon_rad),                  cos(obs_lon_rad),                 0;
        -sin(obs_lat_rad)*cos(obs_lon_rad), -sin(obs_lat_rad)*sin(obs_lon_rad), cos(obs_lat_rad);
         cos(obs_lat_rad)*cos(obs_lon_rad),  cos(obs_lat_rad)*sin(obs_lon_rad), sin(obs_lat_rad)
    ];

    % Line-of-sight vector in ENU
    enu = A_enu * dp_ecef;

    E = enu(1);
    N = enu(2);
    U = enu(3);

    % Ranges
    rho_h = sqrt(E^2 + N^2);
    rho   = sqrt(E^2 + N^2 + U^2);

    % Nonlinear measurement
    az = atan2(E, N);
    az = mod(az, 2*pi);

    el = atan2(U, rho_h);

    z = [
        az;
        el;
        rho;
        v
    ];

    % Position partial derivatives in ECEF
    dp_dr = [
        cos(phi)*cos(lambda);
        cos(phi)*sin(lambda);
        sin(phi)
    ];

    dp_dlambda = [
        -r*cos(phi)*sin(lambda);
         r*cos(phi)*cos(lambda);
         0
    ];

    dp_dphi = [
        -r*sin(phi)*cos(lambda);
        -r*sin(phi)*sin(lambda);
         r*cos(phi)
    ];

    % Convert position partials to ENU partials
    enu_r      = A_enu * dp_dr;
    enu_lambda = A_enu * dp_dlambda;
    enu_phi    = A_enu * dp_dphi;

    Er = enu_r(1);
    Nr = enu_r(2);
    Ur = enu_r(3);

    Elambda = enu_lambda(1);
    Nlambda = enu_lambda(2);
    Ulambda = enu_lambda(3);

    Ephi = enu_phi(1);
    Nphi = enu_phi(2);
    Uphi = enu_phi(3);

    % Avoid singularities directly overhead or at zero range
    eps_val = 1e-12;

    if rho_h < eps_val
        warning('Horizontal range is near zero. Azimuth/elevation Jacobian may be singular.');
        rho_h = eps_val;
    end

    if rho < eps_val
        warning('Slant range is near zero. Range Jacobian may be singular.');
        rho = eps_val;
    end

    az_den = E^2 + N^2;
    if az_den < eps_val
        az_den = eps_val;
    end

    % Azimuth partials
    dAz_dr      = (N*Er      - E*Nr)      / az_den;
    dAz_dlambda = (N*Elambda - E*Nlambda) / az_den;
    dAz_dphi    = (N*Ephi    - E*Nphi)    / az_den;

    % Elevation partials
    dEl_dr = (rho_h*Ur - U*((E*Er + N*Nr)/rho_h)) / rho^2;

    dEl_dlambda = ...
        (rho_h*Ulambda - U*((E*Elambda + N*Nlambda)/rho_h)) / rho^2;

    dEl_dphi = ...
        (rho_h*Uphi - U*((E*Ephi + N*Nphi)/rho_h)) / rho^2;

    % Range partials
    dRho_dr      = (E*Er      + N*Nr      + U*Ur)      / rho;
    dRho_dlambda = (E*Elambda + N*Nlambda + U*Ulambda) / rho;
    dRho_dphi    = (E*Ephi    + N*Nphi    + U*Uphi)    / rho;

    % Linearized measurement matrix
    H = zeros(4,6);

    H(1,:) = [dAz_dr,  dAz_dlambda,  dAz_dphi,  0, 0, 0];
    H(2,:) = [dEl_dr,  dEl_dlambda,  dEl_dphi,  0, 0, 0];
    H(3,:) = [dRho_dr, dRho_dlambda, dRho_dphi, 0, 0, 0];
    H(4,:) = [0,       0,            0,         1, 0, 0];
end

% State: [radius; longitude; latitude; speed; flight_path_angle; heading]
x = [
    6371.8e3 + 30000;   % r, m
    deg2rad(-75);       % longitude, rad
    deg2rad(40);        % latitude, rad
    5000;               % speed, m/s
    deg2rad(-5);        % gamma, rad
    deg2rad(90)         % psi, rad
];

% Observer position
obs_lat_rad = deg2rad(39.95);
obs_lon_rad = deg2rad(-75.16);
obs_alt_m   = 0;

Re = 6371.8e3;
obs_r = Re + obs_alt_m;

obs_ecef_m = [
    obs_r*cos(obs_lat_rad)*cos(obs_lon_rad);
    obs_r*cos(obs_lat_rad)*sin(obs_lon_rad);
    obs_r*sin(obs_lat_rad)
];

[z, H] = hgv_aer_speed_measurement(x, obs_ecef_m, obs_lat_rad, obs_lon_rad);

az_rad    = z(1);
el_rad    = z(2);
range_m   = z(3);
speed_mps = z(4);

az_deg = rad2deg(az_rad);
el_deg = rad2deg(el_rad);

z = [azimuth; elevation; range; speed];

H = dzdx;

x = [r; lambda; phi; v; gamma; psi];