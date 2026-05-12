
function [f_sym, A_sym, B_sym, X, U, P] = hgv_dynamics()
%HGV_DYNAMICS_SYMBOLIC_NO_CORIOLIS
% Symbolic 3-DOF spherical-Earth HGV equations of motion without Coriolis.
%
% State:
%   X = [r; lambda; phi; v; gamma; psi], radius, longitude, latitude,
%   speed, flight path, heading
%
% Control:
%   U = sigma
%
% Parameters:
%   P = [Re; mu; rho0; hs; S; m; CL; CD]

    syms r lambda_sym phi_sym v gamma_sym psi_sym lon_bias_est  lat_bias_est real
    syms sigma_sym real
    syms Re mu rho0 hs S m CL CD D L longitude_bias latitude_bias real 

    X = [
        r;
        lambda_sym;
        phi_sym;
        v;
        gamma_sym;
        psi_sym
        lon_bias_est
        lat_bias_est
    ];

    %U = [sigma_sym; D; L];
     U = [D; L; longitude_bias; latitude_bias];

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
    % D = qbar * CD * S / m;
    % L = qbar * CL * S / m;

    % Equations of motion
    f_sym = sym(zeros(6,1));

    f_sym(1) = v * sin(gamma_sym);

    f_sym(2) = longitude_bias + (v * cos(gamma_sym) * sin(psi_sym)) / (r * cos(phi_sym));

    f_sym(3) = latitude_bias + (v * cos(gamma_sym) * cos(psi_sym)) / r;

    f_sym(4) = -D - g * sin(gamma_sym);

    f_sym(5) = (1/v) * ( ...
                 L * cos(0) ...
               + ((v^2 / r) - g) * cos(gamma_sym));

    f_sym(6) = (1/v) * ( ...
                 (L * sin(0)) / cos(gamma_sym) ...
               + (v^2 / r) * cos(gamma_sym) * sin(psi_sym) * tan(phi_sym));
    f_sym(7) = 0;

    f_sym(8) = 0;


    % Continuous-time linearization
    A_sym = simplify(jacobian(f_sym, X));
    B_sym = simplify(jacobian(f_sym, U));
end