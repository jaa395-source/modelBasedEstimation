function [Hfun_tgt, Hfun_obs] = azElRangeH_symbolic()
% azElRangeH_symbolic
%
% Symbolic linearized measurement matrix for spacecraft 2 measuring
% azimuth, elevation, and range to spacecraft 1.
%
% Measurement:
%   z = [az; el; rho]
%
% Estimated target state:
%   x1 = [r1; phi1; lambda1]
%
% Observer spacecraft 2 is treated as known:
%   x2 = [r2; phi2; lambda2]
%
% Angles are in radians.
%
% Linearized model:
%   dz = H * dx1 + v
%
% Outputs:
%   Hsym = symbolic 3x3 Jacobian d[az, el, rho]/d[r1, phi1, lambda1]
%   Hfun = numeric MATLAB function handle for evaluating Hsym

    syms r1 phi1 lambda1 real
    syms r2 phi2 lambda2 real

    assumeAlso(r1 > 0)
    assumeAlso(r2 > 0)

    % Longitude difference from observer to target
    dlam = lambda1 - lambda2;

    % NED line-of-sight components from spacecraft 2 to spacecraft 1
    N = r1 * ( ...
        cos(phi2)*sin(phi1) ...
        - sin(phi2)*cos(phi1)*cos(dlam) );

    E = r1 * cos(phi1) * sin(dlam);

    D = r2 - r1 * ( ...
        cos(phi2)*cos(phi1)*cos(dlam) ...
        + sin(phi2)*sin(phi1) );

    % Measurement equations
    az  = atan2(E, N);
    el  = atan2(-D, sqrt(N^2 + E^2));
    rho = sqrt(N^2 + E^2 + D^2);

    z = [az; el; rho];

    % Target state
    x1 = [r1; phi1; lambda1];
    x2 = [r2; phi2; lambda2];

    % Linearized H matrix
    Hsym_tgt = simplify(jacobian(z, x1));
    Hsym_obs = simplify(jacobian(z, x2));

    % Numeric evaluator
    Hfun_tgt = matlabFunction(Hsym_tgt, ...
        'Vars', {r1, phi1, lambda1, r2, phi2, lambda2});

    Hfun_obs = matlabFunction(Hsym_obs, ...
        'Vars', {r1, phi1, lambda1, r2, phi2, lambda2});
end