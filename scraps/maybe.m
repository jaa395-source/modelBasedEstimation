function [f_sym, A_sym, A0, f0] = linearize_entry_dynamics(X0, sigma0, params)
% State vector:
% X = [r; theta; phi; v; gamma; psi]
%
% Input:
% sigma = bank angle
%
% params fields:
%   Re, hs, rho0, CD, CL, S, m, mu, Cgamma, Cpsi
%
% h = r - Re
% g = mu / r^2
% rho = rho0 * exp(-h/hs)
% q = rho * v^2

    % Define symbolic variables
    syms r theta phi v gamma psi sigma real
    syms Re hs rho0 CD CL S m mu Cgamma Cpsi real

    % State vector
    X = [r;
         theta;
         phi;
         v;
         gamma;
         psi];

    % Altitude
    h = r - Re;

    % Atmosphere and gravity
    rho = rho0 * exp(-h/hs);
    q = rho * v^2;

    g = mu / r^2;

    % Aerodynamic accelerations
    aD = q * CD * S / m;
    aL = q * CL * S / m;

    % Nonlinear dynamics
    f_sym = sym(zeros(6,1));

    f_sym(1) = v*sin(gamma);

    f_sym(2) = (v*cos(gamma)*sin(psi)) / (r*cos(phi));

    f_sym(3) = (v*cos(gamma)*cos(psi)) / r;

    f_sym(4) = -aD - g*sin(gamma);

    f_sym(5) = (1/v)*( ...
                 aL*cos(sigma) ...
               + ((v^2/r) - g)*cos(gamma) ...
               ) + Cgamma;

    f_sym(6) = (1/v)*( ...
                 (aL*sin(sigma))/cos(gamma) ...
               + (v^2/r)*cos(gamma)*sin(psi)*tan(phi) ...
               ) + Cpsi;

    % Jacobian with respect to state
    A_sym = simplify(jacobian(f_sym, X));

    % Substitute numerical parameter values
    param_syms = [Re, hs, rho0, CD, CL, S, m, mu, Cgamma, Cpsi];

    param_vals = [params.Re, ...
                  params.hs, ...
                  params.rho0, ...
                  params.CD, ...
                  params.CL, ...
                  params.S, ...
                  params.m, ...
                  params.mu, ...
                  params.Cgamma, ...
                  params.Cpsi];

    % Substitute trajectory point X0 and sigma0
    sub_syms = [X;
                sigma;
                param_syms.'];

    sub_vals = [X0(:);
                sigma0;
                param_vals(:)];

    A0 = double(subs(A_sym, sub_syms, sub_vals));
    f0 = double(subs(f_sym, sub_syms, sub_vals));
end

%%
axis_lims = [0 440 -30 30];
plot_estimator2(tvec,xinfo_many,Pinfo_many,x_true,z,[1 4], {"range (m)", "velocity (m/s)"}, axis_lims);
sgtitle('Information Filter (' + msmt_count + ' msmts)','FontWeight','bold','Fontsize',16);

axis_lims = [0 440 -0.7 0.7];
plot_estimator2(tvec,xinfo_many,Pinfo_many,x_true,z,[2 3], {"longitude (rad)", "latitude (rad)"}, axis_lims);
sgtitle('Information Filter (' + msmt_count + ' msmts)','FontWeight','bold','Fontsize',16);

disp("Information filter run for " + msmt_count + " sensors");